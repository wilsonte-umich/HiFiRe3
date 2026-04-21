//! Structure and methods for propagating aggregated alignment 
//! information to files used by the R Shiny app and elsewhere
//! for relating alignments to SV junctions.

// dependencies
use std::error::Error;
use rust_htslib::bam::{record::Record as BamRecord};
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use mdi::{InputCsv, OutputCsv};
use genomex::bam::junction::get_unpacked_node_aln;
use genomex::genome::{Chroms, BedGraphU16};
use crate::sites::SiteMatches;
use super::{JunctionChromWorker, FinalJunction};

// constants
const DESERIALIZATION_ERROR: &str = "Deserialization error when reading alignment segments";

/// An AlignmentSegment describes a single alignment segment within a read.
/// Might be a 5'/proximal, middle, or 3'/distal alignment segment.
/// Unique values are aggregated and counted for tabulation and plotting.
#[derive(PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, Default)]
pub struct AlignmentSegment {
    chrom_index1:  u8,
    ref_pos5:      u32, // 5' site_pos1 for first aln of RE library OR actual 5' ref_pos1 for other alns or non-RE
    ref_pos3:      u32, // actual 3' ref_pos1 of all alignments
    ref_proj3:     u32, // projected 3' site_pos1 for last aln of a RE library OR actual 3' ref_pos1 for other alns or non-RE
    strand_index0: u8,  // 0|1 = +|- = forward|reverse strand
    jxn_types:     u8,  // bit-encoded junction type(s) flanking this alignment segment (one or two)
    n_jxns:        u8,  // number of junctions in the read path containing this alignment segment, including traversal failures
    aln_i:         u8,  // index of this alignment segment within the read path
    n_observed:    u16, // number of times this unique alignment segment was observed
}
impl AlignmentSegment {
    
    /// Create a new AlignmentSegment from a read alignment in a BamRecord.
    pub fn new(
        alns:      &[BamRecord],
        aln_i:     usize,
        jxn_types: u8,
        tool:      &JunctionChromWorker,
    ) -> Self {
        let n_alns = alns.len();
        let max_aln_i = n_alns - 1;

        // initilize site lookup; from_bam_record returns empty SiteMatches if not applicable
        let site_matches = SiteMatches::from_bam_record(&alns[aln_i]);

        // 5' end pos, set to site_pos1 on first alignment in read when expecting outer RE sites
        // -1 correction on reverse aln ref_pos5 accounts for site_pos1 encoding as first base after cleaved bond
        let ref_pos5 = match (tool.expecting_endpoint_re_sites, aln_i){
            (true, 0) => site_matches.site5.pos1 as u32 - if alns[aln_i].is_reverse(){ 1 } else { 0 },
            _ => get_unpacked_node_aln(&alns[aln_i], &tool.header, 5, tool.chroms).2 as u32,
        };

        // actual alignment 3' end included for every alignment in every read path
        let (_chrom, chrom_index1, ref_pos3, is_reverse) = 
            get_unpacked_node_aln(&alns[aln_i], &tool.header, 3, tool.chroms);
        let ref_pos3 = ref_pos3 as u32;

        // projected alignment 3' end, set to site_pos1 on last alignment in read when expecting outer RE sites
        // for single-aln reads, both ref_pos5_observed and ref_pos3_projected might both be RE site_pos1
        // -1 correction on forward aln ref_proj3 accounts for site_pos1 encoding as first base after cleaved bond
        let ref_proj3 = match (tool.expecting_endpoint_re_sites, aln_i){
            (true, i) if i == max_aln_i => site_matches.proj3.pos1 as u32 - if alns[aln_i].is_reverse(){ 0 } else { 1 },
            _ => ref_pos3,
        };

        // return the result
        AlignmentSegment {
            chrom_index1:   chrom_index1 as u8,
            ref_pos5,
            ref_pos3,
            ref_proj3,
            strand_index0:  is_reverse as u8,
            jxn_types,
            n_jxns:         n_alns as u8 - 1,
            aln_i:          aln_i as u8,
            n_observed:     1,
        }
    }

    /// Write a vector of AlignmentSegments to a file, one row per
    /// unique sorted segment with a count.
    pub fn write_sorted(
        mut aln_segments: Vec<AlignmentSegment>, 
        filepath: &str,
        n_cpu:    u32,
    ) -> Result<(), Box<dyn Error>> {
        aln_segments.par_sort_unstable();
        let mut writer = OutputCsv::open_csv(filepath, b'\t', false, Some(n_cpu));
        for chunk in aln_segments.chunk_by_mut(|a, b| a == b) {
            chunk[0].n_observed = chunk.len() as u16;
            writer.serialize(&chunk[0]);
        }
        writer.close();
        Ok(())
    }

    /// Sort and merge distal AlignmentSegments into a set of first_alns files 
    /// partitioned by chromosome and write to a final output file.
    pub fn write_merged(
        mut distal_alns:   Vec<AlignmentSegment>,
        chrom_file_prefix: &str,
        aln_filepath:      &str,
        cvg_filepath:      &str,
        chroms:            &Chroms,
        final_jxns:        &mut Vec<FinalJunction>,
        n_cpu:             u32,
    ) -> Result<(usize, usize, usize, usize), Box<dyn Error>> {

        // initialize iteration over sorted distal alignments
        distal_alns.par_sort_unstable();
        let n_distal_alns = distal_alns.len();
        let mut distal_aln_chunks: Vec<_> = distal_alns.par_chunk_by_mut(|a, b| a == b).collect();
        let n_uniq_distal_alns = distal_aln_chunks.len();
        for distal_aln_chunk in distal_aln_chunks.iter_mut() {
            distal_aln_chunk[0].n_observed = distal_aln_chunk.len() as u16;
        }
        let mut distal_iter = distal_aln_chunks.iter();
        let mut next_distal_aln_chunk = distal_iter.next();

        // intialize output writers and first alignment count
        let mut aln_writer = OutputCsv::open_csv(aln_filepath, b'\t', false, Some(n_cpu));
        let mut cvg_writer = OutputCsv::open_csv(cvg_filepath, b'\t', false, Some(n_cpu));
        let mut n_first_alns: usize = 0;
        let mut n_uniq_first_alns: usize = 0;

        // process each canonical chromosome in order, with previously sorted first alignments from file
        for chrom_index in &chroms.canonical_indices {
            let chrom_index_padded = format!("{:02}", chrom_index);
            let first_alns_file = format!("{}.first_alns.chr{}.gz", chrom_file_prefix, chrom_index_padded);
            if !std::path::Path::new(&first_alns_file).exists() { continue; }
            eprintln!("    {}", chroms.rev_index[chrom_index]);
            let mut first_alns = InputCsv::open_file(&first_alns_file, b'\t', false);

            // initialize coverage map for this chromosome
            let mut coverage_map = vec![0_u16; chroms.index_sizes[&chrom_index] as usize];

            // run the chromosome's first alignments, merge in distal alignments, increment coverage map
            for first_aln in first_alns.deserialize::<AlignmentSegment>() {
                let first_aln = first_aln.expect(DESERIALIZATION_ERROR);
                while let Some(distal_aln_chunk) = next_distal_aln_chunk {
                    let distal_aln = &distal_aln_chunk[0];
                    if distal_aln < &first_aln {
                        if distal_aln.chrom_index1 == *chrom_index { distal_aln.add_to_map(&mut coverage_map); }
                        aln_writer.serialize(distal_aln);
                        next_distal_aln_chunk = distal_iter.next();
                    } else {
                        break;
                    }
                }
                first_aln.add_to_map(&mut coverage_map);
                aln_writer.serialize(&first_aln);
                n_first_alns += first_aln.n_observed as usize;
                n_uniq_first_alns += 1;
            }

            // update final junction breakpoint coverage fields from coverage map
            // known limitation: may miss counting a few distal_alns at very end of chrom
            for jxn in final_jxns.iter_mut(){ 
                if jxn.chrom_index1_1 == *chrom_index {
                    jxn.bkpt_coverage_1 = coverage_map[jxn.ref_pos1_1 as usize - 1];
                }
                if jxn.chrom_index1_2 == *chrom_index {
                    jxn.bkpt_coverage_2 = coverage_map[jxn.ref_pos1_2 as usize - 1];
                }
            }

            // write entries to the begraph coverage file, omitting zero-coverage regions
            let mut chrom_offset0 = 0;
            for coverage_chunk in coverage_map.chunk_by(|a, b| a == b){
                let n_pos = coverage_chunk.len();
                if coverage_chunk[0] > 0 {
                    cvg_writer.serialize(&BedGraphU16{
                        chrom_index: *chrom_index,
                        start0:      chrom_offset0 as u32,
                        end1:        (chrom_offset0 + n_pos) as u32,
                        coverage:    coverage_chunk[0],
                    });
                }
                chrom_offset0 += n_pos;
            }
        }

        // write any remaining distal alignments at the end of the last chromosome
        while let Some(distal_aln_chunk) = next_distal_aln_chunk {
            let distal_aln = &distal_aln_chunk[0];
            aln_writer.serialize(&distal_aln);
            next_distal_aln_chunk = distal_iter.next();
        }

        // close the writer and return the number of first alignments written
        aln_writer.close();
        cvg_writer.close();
        Ok((n_first_alns, n_uniq_first_alns, n_distal_alns, n_uniq_distal_alns))
    }

    /// Increment coverage map for the reference positions spanned by this alignment segment.
    fn add_to_map(&self, coverage_map: &mut [u16]) {
        let (start0, end1) = if self.ref_pos3 > self.ref_pos5 {
            (self.ref_pos5 as usize - 1, self.ref_pos3 as usize)
        } else {
            (self.ref_pos3 as usize - 1, self.ref_pos5 as usize)
        };
        coverage_map[start0..end1].iter_mut()
            .for_each(|pos_cov| *pos_cov = pos_cov.saturating_add(self.n_observed));
    }

    /// Merge two or more alignment segment files.
    pub fn merge_and_write_sorted(
        merge_input_dirs: &[&str],
        filepath: &str,
        ncpu:     u32,
    ) -> Result<(), Box<dyn Error>> {
        let mut aln_segments: Vec<AlignmentSegment> = Vec::new();
        for dir in merge_input_dirs {
            let mut reader = InputCsv::open_file_from_glob(
                dir, "alignments.txt.bgz", 
                b'\t', false
            )?;
            for result in reader.deserialize(){
                let aln_segment = result?;
                aln_segments.push(aln_segment);
            }
        }
        AlignmentSegment::write_sorted(aln_segments, filepath, ncpu)
    }
}

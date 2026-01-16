//! Structure and methods for propagating aggregated alignment information
//! to files used by the R Shiny app and elsewhere.

// dependencies
use std::error::Error;
use std::io::Write;
use rust_htslib::bgzf::Writer as BgzWriter;
use rust_htslib::bam::{record::Record as BamRecord};
use rayon::prelude::*;
use genomex::bam::{tags, junction::get_unpacked_node_aln};
use crate::formats::hf3_tags::*;
use super::JunctionTool;

// constants
const NODE5_OBSERVED_SITE_POS:  usize = 1;
const NODE3_PROJECTED_SITE_POS: usize = 7;

/// An AlignmentSegment describes a single alignment segment within a read.
/// Might be a 5'/proximal, middle, or 3'/distal alignment segment.
/// Unique values are aggregated and counted for tabulation and plotting.
#[derive(PartialEq, Eq, PartialOrd, Ord)]
pub struct AlignmentSegment {
    chrom_index1:       u8,
    ref_pos5_observed:  u32, // 5' site_pos1 for first aln of RE library OR actual 5' ref_pos1 for other alns or non-RE
    ref_pos3_observed:  u32, // actual 3' ref_pos1 of all alignments
    ref_pos3_projected: u32, // (projected) 3' site_pos1 for last aln of a RE library OR actual 3' ref_pos1 for other alns or non-RE
    strand_index0:      u8,  // 0|1 = +|- = forward|reverse strand
    jxn_types:          u8,  // bit-encoded junction type(s) flanking this alignment segment (one or two)
    n_jxns:             u8,  // number of all junctions in the read path containing this alignment segment; includes traversal failures
    aln_i:              u8,  // index of this alignment segment within the read path
}
impl AlignmentSegment {
    
    /// Create a new AlignmentSegment from a read alignment.
    pub fn new(
        alns:      &[BamRecord],
        aln_i:     usize,
        jxn_types: u8,
        tool:      &JunctionTool,
    ) -> Self {
        let n_alns = alns.len();
        let max_aln_i = n_alns - 1;

        // initilize site lookup if applicable
        let closest_sites = if tool.expecting_endpoint_re_sites{
            tags::get_tag_i32_vec(&alns[aln_i], CLOSEST_SITES)
        } else {
            Vec::new()
        };

        // 5' end pos, set to site_pos1 on first alignment in read when expecting outer RE sites
        let ref_pos5_observed = match (tool.expecting_endpoint_re_sites, aln_i){
            (true, 0) => closest_sites[NODE5_OBSERVED_SITE_POS] as u32,
            _ => get_unpacked_node_aln(&alns[aln_i], &tool.header, 5, tool.chroms).2 as u32,
        };

        // actual alignment 3' end included for every alignment in every read path
        let (_chrom, chrom_index1, ref_pos3_observed, is_reverse) = 
            get_unpacked_node_aln(&alns[aln_i], &tool.header, 3, tool.chroms);
        let ref_pos3_observed = ref_pos3_observed as u32;

        // projected alignment 3' end, set to site_pos1 on last alignment in read when expecting outer RE sites
        // for single-aln reads, both ref_pos5_observed and ref_pos3_projected might both be RE site_pos1
        let ref_pos3_projected = match (tool.expecting_endpoint_re_sites, aln_i){
            (true, i) if i == max_aln_i => closest_sites[NODE3_PROJECTED_SITE_POS] as u32,
            _ => ref_pos3_observed,
        };

        // return the result
        AlignmentSegment {
            chrom_index1:   chrom_index1 as u8,
            ref_pos5_observed,
            ref_pos3_observed,
            ref_pos3_projected,
            strand_index0:  is_reverse as u8,
            jxn_types,
            n_jxns:         n_alns as u8 - 1,
            aln_i:          aln_i as u8,
        }
    }

    /// Write a vector of AlignmentSegments to a file, one row per
    /// unique sorted segment with a count.
    pub fn write_sorted_with_count(
        mut aln_segments: Vec<AlignmentSegment>, 
        file_path: &str
    ) -> Result<(), Box<dyn Error>> {
        let mut bgz_writer = BgzWriter::from_path(file_path)?;
        aln_segments.par_sort_unstable();
        for chunk in aln_segments.chunk_by(|a, b| a == b) {
            let seg = &chunk[0];
            writeln!(
                bgz_writer, 
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                seg.chrom_index1,
                seg.ref_pos5_observed,
                seg.ref_pos3_observed,
                seg.ref_pos3_projected,
                seg.strand_index0,
                seg.jxn_types,
                seg.n_jxns,
                seg.aln_i,
                chunk.len()
            )?;
        }
        Ok(())
    }    
}

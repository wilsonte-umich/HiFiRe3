//! Process reads with first alignments on a specific chromosome, 
//! provided as a message on a channel.

// dependencies
use std::error::Error;
use std::fs::read_to_string;
use crossbeam::channel::{Receiver, Sender};
use rust_htslib::bam::{Reader, Read, Record as BamRecord};
use mdi::pub_key_constants;
use mdi::workflow::Config;
use genomex::genome::Chroms;
use genomex::bam::{tags as bam_tags, chrom::get_chrom_index};
use genomex::sam::{junction::JunctionType};
use crate::formats::hf3_tags::*;
use crate::junctions::*;
use crate::inserts::ReadLevelMetadata;

// constants
pub_key_constants!(
    // from environment variables
    EXPECTING_ENDPOINT_RE_SITES
    REJECTING_JUNCTION_RE_SITES
    SEQUENCING_PLATFORM
    INDEX_FILE_PREFIX_WRK
);

// process chromosomes received on the channel
pub fn process_chrom <'a>(
    chroms:    &'a Chroms,
    rx_chrom:  Receiver<(String, u8)>,
    tx_data:   Sender<ChromWorkerData>,
    n_cpu:     u32,
) -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_bool_env(&[EXPECTING_ENDPOINT_RE_SITES, REJECTING_JUNCTION_RE_SITES]);
    cfg.set_string_env(&[SEQUENCING_PLATFORM, INDEX_FILE_PREFIX_WRK]);
    let chrom_file_prefix = cfg.get_string(INDEX_FILE_PREFIX_WRK);

    // process chromosomes received on the channel
    for (chrom_name, chrom_index) in rx_chrom.iter() {

        // open the input BAM file
        // all reads are on-target and have first alignment on chrom
        let chrom_index_padded = format!("{:02}", chrom_index);
        let chrom_bam_path  = format!("{}.chr{}.bam",        chrom_file_prefix, chrom_index_padded);
        let read_count_path = format!("{}.chr{}.read_count", chrom_file_prefix, chrom_index_padded);
        let mut chrom_bam = Reader::from_path(&chrom_bam_path)?;
        let header = chrom_bam.header().clone();
        let chrom_read_count: usize = read_to_string(&read_count_path)?.trim().parse().unwrap();

        // guesstimate the expected number of alignments for pre-allocation
        let expecting_endpoint_re_sites = *cfg.get_bool( EXPECTING_ENDPOINT_RE_SITES);
        let first_alns_capacity = if expecting_endpoint_re_sites {
            chrom_read_count / 4 // RE-digested inserts are expected to have fewer unique first alignments
        } else {
            chrom_read_count     // uses 19 bytes per read, so 100M reads requires 1.9GB RAM, more than expected per chrom
        };

        // assemble the chromosome worker tool
        let mut tool = JunctionChromWorker {
            expecting_endpoint_re_sites,
            chroms,
            header,
            first_alns: Vec::with_capacity(first_alns_capacity), // pre-allocate for first alignments
            tx_data:    &tx_data,
        };

        // short-circuit empty chromosomes with no alignments
        let mut records = chrom_bam.records();
        let first_aln = records.next();
        if first_aln.is_none() {
            tx_data.send(ChromWorkerData::ChromReadCount((chrom_name, 0)))?;
            continue;
        }

        // process alignment records by read
        // all alignments are in query read order
        let mut alns: Vec<BamRecord> = vec![first_aln.unwrap()?];
        let mut chrom_read_count:usize = 0;
        for result in records {
            let record = result?;
            if record.qname() != alns[0].qname() { 
                process_alns(&alns, &mut tool)?;
                chrom_read_count += 1;
                alns.clear();
            }
            alns.push(record);
        }
        chrom_read_count += 1;
        process_alns(&alns, &mut tool)?;

        // sort and print first alignments per chromosome to limit memory overhead (merged later)
        let first_alns_file = format!("{}.first_alns.chr{}.gz", chrom_file_prefix, chrom_index_padded);
        AlignmentSegment::write_sorted(
            tool.first_alns, 
            &first_alns_file,
            n_cpu,
        )?;

        // send chrom read count to main thread
        tx_data.send(ChromWorkerData::ChromReadCount((chrom_name, chrom_read_count)))?;
    }
    Ok(())
}

// process all alignments for a single read
fn process_alns(
    alns:   &[BamRecord], 
    tool:   &mut JunctionChromWorker,
) -> Result<(), Box<dyn Error>>{

    // skip reads with any alignment to a non-canonical chromosome
    if alns.iter().any(|aln| {
        get_chrom_index(aln, &tool.header, &tool.chroms, false).is_none()
    }) { return Ok(()); }

    // get read-level metadata from first alignment
    let n_alns = alns.len();
    let max_aln_i = n_alns - 1;
    let sample_bit = bam_tags::get_tag_u32(&alns[0], SAMPLE_BIT);
    let read_data = ReadLevelMetadata::from_bam_records(alns, sample_bit);
    let qlen = alns[0].seq_len() as u32;

    // has_committed_jxn is not the same as READ_HAS_PASSED_JXN or n_alns > 1
    // traversal failures, and only traversal failures, are skipped here; chimeras persist
    let read_has_jxn = n_alns > 1;
    let mut has_committed_jxn = false; 
    let mut prev_jxn_type: u8 = JunctionType::Proper as u8;

    // run the read's alns
    for aln_i in 0..=max_aln_i {
        let mut this_jxn_type: u8 = JunctionType::Proper as u8;

        // assess any junctions
        if read_has_jxn && aln_i < max_aln_i {
            let (ordered_jxn, jxn_instance) = extract_junction(
                aln_i as u8,
                &alns[aln_i],
                &alns[aln_i + 1],
                &read_data,
                qlen,
            );
            this_jxn_type = ordered_jxn.jxn_type;
            let jxn_failure_flag = bam_tags::get_tag_u8_default(&alns[aln_i], JXN_FAILURE_FLAG, 0);
            let passed_traversal_delta = jxn_failure_flag & JxnFailureFlag::TraversalDelta as u8 == 0;
            if passed_traversal_delta { // only commit junctions that passed traversal delta
                tool.tx_data.send(ChromWorkerData::Junction((ordered_jxn, jxn_instance)))?;
                has_committed_jxn = true;
            }
        }

        // send assembled alignment segments to first_alns Vec or distal_alns main thread channel
        let aln_segment = AlignmentSegment::new(
            alns,
            aln_i,
            prev_jxn_type | this_jxn_type, // jxn_types flanking this aln (could be two types for a middle aln)
            tool,
        );
        if aln_i == 0 {
            tool.first_alns.push(aln_segment);
        } else {
            tool.tx_data.send(ChromWorkerData::DistalAln(aln_segment))?;
        }
        prev_jxn_type = this_jxn_type;
    }

    // send SV read paths to main thread channel
    if has_committed_jxn {
        let read_path = SvReadPath::from_alns(alns, &tool.header);
        tool.tx_data.send(ChromWorkerData::ReadPath(read_path))?;
    }
    Ok(())
}

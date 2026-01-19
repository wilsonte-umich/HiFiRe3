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
use genomex::bam::tags;
use genomex::sam::{junction::JunctionType};
use crate::formats::hf3_tags::*;
use crate::junctions::*;
use crate::inserts::ReadLevelMetadata;

// constants
pub_key_constants!(
    // from environment variables
    INDEX_FILE_PREFIX_WRK
    SEQUENCING_PLATFORM
    EXPECTING_ENDPOINT_RE_SITES
    REJECTING_JUNCTION_RE_SITES
    // derived variables
    IS_ONT
    // counter keys
    N_READS
    N_READS_ON_TARGET
    N_UNIQ_ALNS
    N_ALNS
    N_REF_BASES
    N_READ_BASES
    N_READS_BY_GENOME
    N_REF_BASES_BY_GENOME
    N_READ_BASES_BY_GENOME
);

// process chromosomes received on the channel
pub fn process_chrom <'a>(
    chroms:        &'a Chroms,
    rx_chrom:      Receiver<(String, u8)>,
    tx_jxn:        Sender<(OrderedJunction, JunctionInstance)>,
    tx_read_path:  Sender<SvReadPath>,
    tx_distal_aln: Sender<AlignmentSegment>,
) -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[INDEX_FILE_PREFIX_WRK, SEQUENCING_PLATFORM]);
    cfg.set_bool_env(&[EXPECTING_ENDPOINT_RE_SITES, REJECTING_JUNCTION_RE_SITES]);
    let chrom_file_prefix = cfg.get_string(INDEX_FILE_PREFIX_WRK);

    // process chromosomes received on the channel
    for (_chrom_name, chrom_index) in rx_chrom.iter() {

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
            chrom_read_count     // uses 17 bytes per read, so 100M reads requires 1.7GB RAM, more than expected per chrom
        };

        // assemble the tool
        let mut tool = JunctionTool {
            expecting_endpoint_re_sites,
            chroms,
            header,
            first_alns:    Vec::with_capacity(first_alns_capacity), // pre-allocate for first alignments
            tx_jxn:        &tx_jxn,
            tx_read_path:  &tx_read_path,
            tx_distal_aln: &tx_distal_aln,
        };

        // process alignment records by read
        // all alignments are in query read order
        let mut records = chrom_bam.records();
        let mut alns: Vec<BamRecord> = vec![records.next().unwrap()?]; // expect there to always be data
        for result in records {
            let record = result?;
            if record.qname() != alns[0].qname() { 
                process_alns(&alns, &mut tool)?;
                alns.clear();
            }
            alns.push(record);
        }
        process_alns(&alns, &mut tool)?;

        // sort and print first alignments per chromosome to limit memory overhead (merged later)
        let first_alns_file = format!("{}.chr{}.first_alns.bgz", chrom_file_prefix, chrom_index_padded);
        AlignmentSegment::write_sorted_with_count(
            tool.first_alns, 
            &first_alns_file
        )?;

    }
    Ok(())
}

// process all alignments for a single read
fn process_alns(
    alns:   &[BamRecord], 
    tool:   &mut JunctionTool,
) -> Result<(), Box<dyn Error>>{
    let n_alns = alns.len();
    let max_aln_i = n_alns - 1;

    // get read-level metadata from first alignment
    let read_data = ReadLevelMetadata::from_bam_records(alns);

    // has_committed_jxn is not the same as READ_HAS_PASSED_JXN or n_alns > 1
    // traversal failures, and only traversal failures, are skipped here; chimeras persist
    let read_has_jxn = n_alns > 1;
    let mut has_committed_jxn = false; 
    let mut prev_jxn_type: u8 = JunctionType::Proper as u8;

    // run the alns
    for aln_i in 0..=max_aln_i {
        let mut this_jxn_type: u8 = JunctionType::Proper as u8;

        // handle the junctions
        if read_has_jxn && aln_i < max_aln_i {
            let (ordered_jxn, jxn_instance) = get_junction(
                &alns[aln_i],
                &alns[aln_i + 1],
                &read_data
            );
            this_jxn_type = ordered_jxn.jxn_type;
            let jxn_failure_flag = tags::get_tag_u8_default(&alns[aln_i], JXN_FAILURE_FLAG, 0);
            let passed_traversal_delta = jxn_failure_flag & JxnFailureFlag::TraversalDelta as u8 == 0;
            if passed_traversal_delta { // only record juctions that passed traversal delta
                tool.tx_jxn.send((ordered_jxn, jxn_instance))?;
                has_committed_jxn = true;
            }
        }

        // handle the alignment segments
        let aln_segment = AlignmentSegment::new(
            alns,
            aln_i,
            prev_jxn_type | this_jxn_type, // jxn_types flanking this aln (could be two types for a middle aln)
            tool,
        );
        if aln_i == 0 {
            tool.first_alns.push(aln_segment);
        } else {
            tool.tx_distal_aln.send(aln_segment)?;
        }
        prev_jxn_type = this_jxn_type;
    }

    // print SV read paths
    if has_committed_jxn {
        let read_path = SvReadPath::from_alns(alns, &tool.header, tool.chroms);
        tool.tx_read_path.send(read_path)?;
    }
    Ok(())
}

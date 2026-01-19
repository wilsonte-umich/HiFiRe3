//! Support for insert comparisons, including deduplication.

// dependencies
use std::str::from_utf8_unchecked;
use rust_htslib::bam::Record as BamRecord;
use genomex::bam::tags;
use genomex::sam::SamRecord;
use crate::formats::hf3_tags::{OUTER_NODES, CHANNEL, INSERT_SIZE};

/// A UniqueInsertSpan holds (re)ordered insert outer nodes,
/// and ONT channel information for deduplication.
/// 
/// Thus, an instance identifies "the canonical orientation of 
/// two outer nodes (for an insert on a specific ONT channel)", 
/// given that inserts/junctions detected in different channels 
/// are considered distinct for deduplication.
/// 
/// UniqueInsertSpan is a read, not an alignment-level
/// property and only needs to be assessed once per read.
#[derive(Hash, Eq, PartialEq, Clone, Copy)]
pub struct UniqueInsertSpan {
    pub node1:   isize,
    pub node2:   isize,
    pub channel: u32,
}
impl UniqueInsertSpan {
    /// Get a unique insert span from a BamRecord as (UniqueInsertSpan, was_reordered), 
    /// where `was_reordered` indicates if the nodes were reordered to achieve 
    /// the canonical representation.
    pub fn from_bam_record(aln: &BamRecord) -> (UniqueInsertSpan, bool) {
        let ordered_outer_nodes = SamRecord::sam_tag_to_paired_nodes(
            &tags::get_tag_str(aln, OUTER_NODES), 
            true
        );
        let channel = tags::get_tag_u32_default(aln, CHANNEL, 0);
        (
            UniqueInsertSpan {
                node1: ordered_outer_nodes.node1,
                node2: ordered_outer_nodes.node2,
                channel,
            },
            ordered_outer_nodes.was_reordered
        )
    }
}

/// ReadLevelMetadata holds metadata that is the same for all alignments of a read.
#[derive(Clone)]
pub struct ReadLevelMetadata {
    pub qname:       String,
    pub insert_size: i32, // INSERT_SIZE tag value
    pub insert_span: UniqueInsertSpan, // reoriented outer nodes and channel for deduplication; reflect site_pos1 when applicable
    pub n_jxns:      u8,  // number of junctions in the read path containing this alignment segment
}
impl ReadLevelMetadata {
    /// Create ReadLevelMetadata from the first BamRecord of a set of BamRecords.
    pub fn from_bam_records(alns: &[BamRecord]) -> Self {
        let (insert_span, _was_reordered) = UniqueInsertSpan::from_bam_record(&alns[0]);
        ReadLevelMetadata {
            qname:        unsafe{ from_utf8_unchecked(alns[0].qname()).to_string() },
            insert_size:  tags::get_tag_i32_default(&alns[0], INSERT_SIZE, alns[0].seq_len() as i32),
            insert_span:  insert_span,
            n_jxns:       alns.len() as u8 - 1,
        }
    }
}

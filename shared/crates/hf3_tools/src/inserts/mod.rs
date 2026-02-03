//! Support for insert comparisons, including deduplication.

// dependencies
use std::str::from_utf8_unchecked;
use rust_htslib::bam::Record as BamRecord;
use genomex::bam::tags;
use genomex::sam::SamRecord;
use crate::formats::hf3_tags::{OUTER_NODES, CHANNEL, INSERT_SIZE};

/// ReadFailureFlag holds bit-encoded data indicating why a read
/// was rejected as a failed read and not used for variant analysis.
#[repr(u8)]
pub enum ReadFailureFlag {
    Unmapped      = 1, 
    OffTarget     = 2,
    SiteLookup    = 4,
    ClipTolerance = 8,
    SiteDistance  = 16,
}

/// A UniqueInsertSpan holds insert outer nodes in the canonical order,
/// and ONT channel information for deduplication.
/// 
/// Thus, an instance identifies "the canonical orientation of two outer 
/// nodes (for an insert on a specific ONT channel)", given that the 
/// same inserts detected in different channels are considered distinct 
/// for deduplication.
/// 
/// UniqueInsertSpan is a read not an alignment-level property and only 
/// needs to be assessed once per read.
#[derive(Hash, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
pub struct UniqueInsertSpan {
    pub node1:   isize, // (re)ordered nodes in the canonical orientation
    pub node2:   isize,
    pub channel: u32,   // ONT channel number from CHANNEL tag; 0 if not applicable
}
impl UniqueInsertSpan {
    /// Get a unique insert span from a BamRecord outer nodes tag as 
    /// (UniqueInsertSpan, was_reordered), 
    /// where `was_reordered` indicates if the nodes were reordered to achieve 
    /// the canonical orientation.
    pub fn from_bam_record(aln: &BamRecord) -> (UniqueInsertSpan, bool) {
        let ordered_outer_nodes = SamRecord::sam_tag_to_paired_nodes(
            &tags::get_tag_str(aln, OUTER_NODES), 
            true
        );
        (
            UniqueInsertSpan {
                node1:   ordered_outer_nodes.node1,
                node2:   ordered_outer_nodes.node2,
                channel: tags::get_tag_u32_default(aln, CHANNEL, 0),
            },
            ordered_outer_nodes.was_reordered
        )
    }
}

/// A ChannelAlignment holds (re)ordered alignment indices, (re)ordered 
/// outer nodes, and ONT channel information for deduplication.
/// 
/// Thus, an instance identifies "the nth alignment of a read after 
/// (re)ordering into the canonical orientation of the two outer nodes 
/// (for an insert on a specific ONT channel)".
#[derive(Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct ChannelAlignment {
    pub unique_insert_span: UniqueInsertSpan,
    pub aln_i:              u8, // up to 256 alignments per read
}

/// ReadLevelMetadata holds metadata that is the same for all alignments and
/// all junctions of a read.
/// 
/// ReadLevelMetadata is a set of read not an alignment-level properties and  
/// only needs to be assessed once per read.
#[derive(Clone)]
pub struct ReadLevelMetadata {
    pub qname:        String, // BAM query name
    pub insert_size:  i32,    // INSERT_SIZE tag value
    pub node1:        isize,  // (re)ordered nodes in the canonical orientation
    pub node2:        isize,
    pub channel:      u32,    // ONT channel number from CHANNEL tag; 0 if not applicable
    pub n_jxns:       u8,     // number of junctions in the read path
    pub is_duplicate: bool,   // whether the read is marked as a duplicate in the BAM flag
}
impl ReadLevelMetadata {
    /// Create ReadLevelMetadata from the first BamRecord of a set of BamRecords.
    pub fn from_bam_records(alns: &[BamRecord]) -> Self {
        let ordered_outer_nodes = SamRecord::sam_tag_to_paired_nodes(
            &tags::get_tag_str(&alns[0], OUTER_NODES), 
            true
        );
        ReadLevelMetadata {
            qname:        unsafe{ from_utf8_unchecked(alns[0].qname()).to_string() },
            insert_size:  tags::get_tag_i32_default(&alns[0], INSERT_SIZE, alns[0].seq_len() as i32),
            node1:        ordered_outer_nodes.node1,
            node2:        ordered_outer_nodes.node2,
            channel:      tags::get_tag_u32_default(&alns[0], CHANNEL, 0),
            n_jxns:       alns.len() as u8 - 1,
            is_duplicate: false, // set later, duplicate status not yet present in BAM flag
        }
    }
}

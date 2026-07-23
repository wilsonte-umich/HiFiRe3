//! Support for collecting read data organized by their source RE fragments.

// dependencies
use std::str::from_utf8_unchecked;
use rustc_hash::FxHashMap;
use rust_htslib::bam::Record as BamRecord;
use genomex::sequence::rc_acgtn_str;
use genomex::bam::tags as bam_tags;
use crate::formats::hf3_tags::*;
use super::*;

/// A unique read span on a known chromosome corresponding to a RE fragment. For 
/// RE-based PacBio sequencing only a relatively limited number of unique read 
/// spans are expected.
#[derive(PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReFragment {
    pub start0: ChromPos0, // BED half-open coordinates (not site positions)
    pub end1:   ChromPos1,
}
impl ReFragment {

    /// Return the standardized ReFragment that matches a newly encountered read.
    pub fn new( 
        mut site5_pos1: ChromPos1, 
        mut site3_pos1: ChromPos1,
        is_reverse: bool,
    ) -> Self {
        if is_reverse { // reverse endpoints for rare reverse strand reads
            (site5_pos1, site3_pos1) = (site3_pos1, site5_pos1);
        }
        ReFragment{
            start0: site5_pos1 - 1, // converts 1-based to 0-based position
            end1:   site3_pos1 - 1, // converts site pos1 AFTER the fragment to fragment end1
        }
    }
}

/// A ReadInstance collects only the essential information from each encountered 
/// read for later re-analysis.
pub struct ReadInstance {
    pub sample_bit: SampleBit,
    pub qname:      String,
    pub seq_bytes:  Vec<u8>, // htslib doesn't make it easy to cache as encoded seq
    pub is_reverse: bool,
    pub cs:     String,
    pub start0: u32, // where cs aln starts and stop on ref, not necessarily at site pos
    pub end1:   u32,
}
impl ReadInstance {

    /// Create a new ReadInstance from its BamRecord.
    pub fn new(aln: &BamRecord) -> Self {
        ReadInstance {
            sample_bit: bam_tags::get_tag_u32(aln, SAMPLE_BIT),
            qname: unsafe { from_utf8_unchecked(aln.qname()).to_string() },
            seq_bytes: aln.seq().as_bytes(),
            is_reverse: aln.is_reverse(),
            cs: bam_tags::get_tag_str(aln, DIFFERENCE_STRING),
            start0: aln.pos() as u32,
            end1: aln.cigar().end_pos() as u32,
        }
    }

    /// Return the ACGTN base sequence of a read on the top reference strand.
    pub fn get_top_strand_seq(&self) -> String{
        let mut seq = unsafe { 
            from_utf8_unchecked(&self.seq_bytes).to_string() 
        };        
        if self.is_reverse { // reverse sequence for rare reverse strand reads
            seq = rc_acgtn_str(&seq);
        }
        seq
    }
}

/// FragmentReads collects encountered ReadInstances for a given ReFragment.
pub struct FragmentReads{
    pub instances: FxHashMap<ReFragment, Vec<ReadInstance>>
}
impl FragmentReads {

    /// Create a new FragmentReads HashMap.
    pub fn new() -> Self{
        Self{
            instances: FxHashMap::default()
        }
    }
    
    /// Add a ReadInstance to the FragmentReads HashMap.
    pub fn insert(&mut self, re_fragment: ReFragment, read_instance: ReadInstance) {
        self.instances
            .entry(re_fragment)
            .or_insert_with(Vec::new)
            .push(read_instance);          
    }
}

/// ReadMap collects information of the specific ReadInstances that reported a 
/// given ReferenceVariant.
pub struct ReadMap {
    pub n_matching_reads:  usize,
    pub read_map: Vec<bool>,
}
impl ReadMap {
    /// Create a new ReadMap.
    pub fn new(n_reads: usize) -> Self{
        Self{
            n_matching_reads: 0,
            read_map: Vec::with_capacity(n_reads)
        }
    } 
}

/// FragmentReferenceVariants collects encountered ReferenceVariants for a given 
/// ReFragment over all FragmentReads.
pub struct FragmentReferenceVariants {
    pub n_reads: usize,
    pub variant_map: FxHashMap<ReferenceVariant, ReadMap>,
}
impl FragmentReferenceVariants {

    /// Create a new FragmentReferenceVariants map.
    pub fn new() -> Self{
        Self{
            n_reads: 0,
            variant_map: FxHashMap::default()
        }
    }

    /// Add one ReferenceVariant from a read to the FragmentReferenceVariants map.
    pub fn insert(
        &mut self, 
        reference_variant: ReferenceVariant,
        read_i: ReadIndex,
    ) {
        let map = self.variant_map
            .entry(reference_variant)
            .or_insert_with(|| ReadMap::new(self.n_reads));
        map.n_matching_reads += 1;
        map.read_map[read_i] = true;        
    }
}

/// FragmentHaplotypeVariants collects encountered HaplotypeVariants for a given 
/// ReFragment over all FragmentReads matching a given haplotype.
pub struct FragmentHaplotypeVariants {
    pub n_reads: usize,
    pub variant_map: FxHashMap<HaplotypeVariant, ReadMap>,
}
impl FragmentHaplotypeVariants {

    /// Create a new FragmentHaplotypeVariants map.
    pub fn new() -> Self{
        Self{
            n_reads: 0,
            variant_map: FxHashMap::default()
        }
    }

    /// Add one HaplotypeVariant from a read to the FragmentHaplotypeVariants map.
    pub fn insert(
        &mut self, 
        haplotype_variant: HaplotypeVariant,
        read_i: ReadIndex,
    ) {
        let map = self.variant_map
            .entry(haplotype_variant)
            .or_insert_with(|| ReadMap::new(self.n_reads));
        map.n_matching_reads += 1;
        map.read_map[read_i] = true;        
    }
}

//! Support for collecting read data organized by their source RE fragments.

// imports
use std::str::from_utf8_unchecked;
use rustc_hash::FxHashMap;
use rust_htslib::bam::Record as BamRecord;
use serde::Serialize;
use genomex::bam::tags as bam_tags;
use crate::formats::hf3_tags::*;
use super::*;

/// A unique read span on a known chromosome corresponding to a RE fragment. For 
/// RE-based PacBio sequencing only a relatively limited number of unique read 
/// spans are expected.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
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
            start0: site5_pos1 - 1, // converts 1-based to 0-based RE fragment start
            end1:   site3_pos1 - 1, // converts site pos1 AFTER the fragment to fragment end1
        }
    }
}

/// A ReadInstance collects only the essential information from each encountered 
/// read for later re-analysis.
pub struct ReadInstance {
    pub sample_bit: SampleBit,
    pub qname:      QName,
    pub seq_bytes:  Vec<BaseByteACGTN>, // htslib doesn't make it easy to cache as encoded seq
    pub qual_bytes: Vec<PhredQual>,
    pub is_reverse: bool,
    pub cs: String,
    pub dd: String,
    pub ref_pos0: ChromPos0, // where cs aln starts on ref, not necessarily at site pos
    pub qry_pos0: SeqPos0,   // where cs aln starts on read AFTER re-orientation to ref strand
}
impl ReadInstance {

    /// Create a new ReadInstance from its BamRecord.
    pub fn new(aln: &BamRecord) -> Self {
        let is_reverse = aln.is_reverse();
        let qry_pos0 = if is_reverse {
            aln.cigar().trailing_softclips()
        } else {
            aln.cigar().leading_softclips()
        } as u32;
        ReadInstance {
            sample_bit: bam_tags::get_tag_u32(aln, SAMPLE_BIT),
            qname: unsafe { from_utf8_unchecked(aln.qname()).to_string() },
            seq_bytes:  aln.seq().as_bytes(), // SEQ and QUAL are not yet re-oriented if is_reverse (done by get_top_strand)
            qual_bytes: aln.qual().to_vec(),  // defer these slower ops that won't be needed on discarded reads
            is_reverse,
            cs: bam_tags::get_tag_str(aln, DIFFERENCE_STRING),  // cs tag is implicity oriented to top strand
            dd: bam_tags::get_tag_str(aln, STRAND_DIFFERENCES), // dd tag is NOT yet re-oriented (done by get_dd_mask)
            ref_pos0: aln.pos() as u32, // leftmost reference position is implicity oriented to top strand
            qry_pos0, // qry_pos0 is already re-oriented if is_reverse (done above)
        }
    }

    /// Return the QUAL of a read on the top reference strand. 
    pub fn get_top_strand_qual(&self) -> Vec<PhredQual> {
        let mut qual = self.qual_bytes.to_vec();
        if self.is_reverse { qual.reverse(); }
        qual
    }
}

/// FragmentReads collects encountered ReadInstances for a given ReFragment.
/// One FragmentReads object is instantiated by each SnvChromWorker.
pub struct FragmentReads{
    pub instances: FxHashMap<ReFragment, Vec<ReadInstance>>
}
impl FragmentReads {

    /// Create a new FragmentReads HashMap.
    pub fn new() -> Self{
        let mut instances = FxHashMap::default();
        instances.reserve(4096); // about 16 Mb of ReFragments to start out
        Self{instances}
    }
    
    /// Add a ReadInstance to the FragmentReads HashMap.
    pub fn insert(&mut self, re_fragment: ReFragment, read_instance: ReadInstance) {
        self.instances
            .entry(re_fragment)
            .or_insert_with(|| Vec::with_capacity(8))
            .push(read_instance);          
    }
}

/// A ReadMapEntry carries bits of information a specific Variant in a specific
/// ReadInstance, including whether it was observed there and at what quality.
#[derive(Clone, Copy)]
pub struct ReadMapEntry {
    has_var:  bool, // immutable record of whether a read reported a specific variant
    pub hap0: bool, // mutable haplotype initally == `has_var`, but may be flipped
    pub avg_qual: PhredQual,
}
impl ReadMapEntry{
    /// Create a new empty ReadMapEntry.
    pub fn new() -> Self{
        Self { 
            has_var:  false, 
            hap0:     false,
            avg_qual: 0 
        }
    }
    /// Get the immutable `has_var` value of a ReadMapEntry.
    pub fn has_var(&self) -> bool {
        self.has_var
    }
}

/// ReadMap collects information of the specific ReadInstances that reported a 
/// given Variant. Allocation is one-time fixed.
pub struct ReadMap {
    pub n_matching_reads:  usize,
    pub read_map: Vec<ReadMapEntry>,
}
impl ReadMap {
    /// Create a new ReadMap.
    pub fn new(n_reads: usize) -> Self{
        Self{
            n_matching_reads: 0,
            read_map: vec![ReadMapEntry::new(); n_reads],
        }
    } 
}

/// FragmentVariants collects encountered ReferenceVariants for a given 
/// ReFragment over all FragmentReads. One FragmentVariants object is 
/// instantiated per SnvChromWorker that is reset as needed per ReFragment.
pub struct FragmentVariants {
    pub n_reads: usize,
    pub variant_map: FxHashMap<Variant, ReadMap>,
}
impl FragmentVariants {

    /// Create a new FragmentVariants map.
    pub fn new() -> Self{
        let mut variant_map = FxHashMap::default();
        variant_map.reserve(128);
        Self{
            n_reads: 0,
            variant_map,
        }
    }

    /// Reset a FragmentVariants map to initialize collection of a new set of 
    /// ReFragment variants.
    pub fn reset(&mut self, n_reads: usize){
        self.n_reads = n_reads;
        self.variant_map.clear();
    }

    /// Add one Variant from a read to the FragmentVariants map.
    pub fn insert(
        &mut self, 
        variant:  Variant,
        read_i:   ReadIndex,
        avg_qual: PhredQual,
    ) {
        let map = self.variant_map
            .entry(variant.clone())
            .or_insert_with(|| ReadMap::new(self.n_reads));
        map.n_matching_reads += 1;
        map.read_map[read_i] = ReadMapEntry{
            has_var: true,
            hap0:    true,
            avg_qual,
        };
    }
}

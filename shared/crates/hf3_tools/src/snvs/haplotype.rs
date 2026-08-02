//! Support for tracking haplotypes and consensuses during SNV calling.

// imports
use serde::{Serialize, Serializer};
use rustc_hash::FxHashMap;
use super::*;

/// Haplotype enumerates the four allowed haplotype values in records.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Debug)]
#[repr(u8)]
pub enum Haplotype {
    Unspecified = 0, // used before haplotypes are resolved
    Haplotype1  = 1,
    Haplotype2  = 2,
    Homozygous  = 3, // thus, bit encoded as 1 + 2
}
/// Helper function to serialize Haplotype as u8.
pub fn serialize_haplotype<S: Serializer>(
    h: &Haplotype, 
    serializer: S
) -> Result<S::Ok, S::Error>{
    serializer.serialize_u8(*h as u8)
}

/// Cache haplotype consensus sequences prior to final printing of read
/// encodings. Also used to store reference sequences prior to consensus build.
pub struct HaplotypeConsensuses {
    pub cache: FxHashMap<(ReFragment, Haplotype), (UppercaseACGTN, Option<String>)>
}
impl HaplotypeConsensuses {

    /// Create a new HaplotypeConsensuses. One cache is created per 
    /// SnvChromWorker.
    pub fn new() -> Self {
        let mut cache = FxHashMap::default();
        cache.reserve(8192);
        Self{ cache }
    }

    /// Add a haplotype consensus or reference sequence to the print cache.
    pub fn insert(
        &mut self,
        re_fragment: &ReFragment,
        haplotype:   Haplotype,
        seq:         UppercaseACGTN, 
        hap_vs_ref:  Option<String>, // None when recording a reference sequence
    ){
        self.cache.insert( 
            (*re_fragment, haplotype), 
            (seq, hap_vs_ref) 
        );
    }

    /// Add a reference sequence to the print cache.
    pub fn insert_reference(
        &mut self,
        tool:        &SnvAnalysisTool,
        worker:      &SnvChromWorker,
        re_fragment: &ReFragment,
    ){
        let seq = tool.fa.view(
            worker.chrom_tid, 
            re_fragment.start0 as usize, 
             re_fragment.end1 as usize
        ).ok()
            .map(|v| v.to_string().to_ascii_uppercase())
            .unwrap_or_else(|| "NA".to_string());
        self.insert(
            re_fragment, 
            Haplotype::Unspecified, 
            seq, 
            None
        );
    }

    /// Get a (sequence, hap_vs_ref) tuple from the cache.
    pub fn get(
        &self, 
        re_fragment: &ReFragment,
        haplotype:   Haplotype,
    ) -> &(UppercaseACGTN, Option<String>) {
        self.cache
            .get(&(*re_fragment, haplotype))
            .expect("Failed to get sequence from HaplotypeConsensuses cache.")
    }
}

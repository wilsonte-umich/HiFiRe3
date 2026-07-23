/// Support for calling and counting specific variants from error-corrected reads.

// dependencies
use serde::Serialize;

/// ReferenceVariant encodes a specific SNV or indel, or a series of operations,
/// observed beginning at a single reference position on a known chromosome.
/// 
/// A ReferenceVariant allows any number of reference bases to be replaced by  
/// any number of non-reference bases, so it is equally capable of representing 
/// substitutions, insertions, deletions, and complex indels. 
/// 
/// ReferenceVariants are an interim stage, eventually replaced if persistent by
/// HaplotypeVariants called relative to a haplotype consensus.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub struct ReferenceVariant {
    // leftmost position when n_ref_bases > 0, or the position immediately preceding an insertion
    pub ref_pos0:    u32,
    // for substitutions and deletions, the number of reference bases replaced by alt_bases         
    pub n_ref_bases: u32,
    // for substitutions and insertions, the non-reference bases replacing the reference bases    
    pub alt_bases:   Option<String>,
}
impl ReferenceVariant {
    /// Create a new ReferenceVariant instance with the specified fields.
    pub fn new(ref_pos0: u32, n_ref_bases: u32, alt_bases: &str) -> Self {
        ReferenceVariant {
            ref_pos0,
            n_ref_bases,
            alt_bases: if alt_bases.is_empty() {
                None 
            } else { 
                Some(alt_bases.to_string()) 
            },
        }
    }

    // /// Get the signed difference in ref vs. alt length.
    // pub fn alt_minus_ref(&self) -> i32 {
    //     self.alt_bases.as_ref().map_or(0, |alt| alt.len() as i32) - 
    //     self.n_ref_bases as i32
    // }
}

/// HaplotypeVariant encodes a specific SNV or indel, or a series of operations,
/// observed beginning at a single position on a resolved fragment haplotype.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub struct HaplotypeVariant {
    pub haplotype: u8, // 1 or 2
    // leftmost position when n_hap_bases > 0, or the position immediately preceding an insertion
    pub hap_pos0:    u32,
    // for substitutions and deletions, the number of haplotype bases replaced by alt_bases         
    pub n_hap_bases: u32,
    // for substitutions and insertions, the non-haplotype bases replacing the haplotype bases    
    pub alt_bases:   Option<String>,
}
impl HaplotypeVariant {
    /// Create a new HaplotypeVariant instance with the specified fields.
    pub fn new(haplotype: u8, hap_pos0: u32, n_hap_bases: u32, alt_bases: &str) -> Self {
        HaplotypeVariant {
            haplotype,
            hap_pos0,
            n_hap_bases,
            alt_bases: if alt_bases.is_empty() {
                None 
            } else { 
                Some(alt_bases.to_string()) 
            },
        }
    }

    // /// Get the signed difference in ref vs. alt length.
    // pub fn alt_minus_ref(&self) -> i32 {
    //     self.alt_bases.as_ref().map_or(0, |alt| alt.len() as i32) - 
    //     self.n_ref_bases as i32
    // }
}

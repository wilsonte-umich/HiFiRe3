/// Support for calling and counting specific variants from error-corrected reads.

// imports
use serde::Serialize;
use super::*;

/// A Variant encodes a specific SNV or indel, or a series of operations,
/// observed beginning at a single reference position on a known chromosome
/// or on a single resolved haplotype.
/// 
/// A Variant allows any number of reference bases to be replaced by any number 
/// of non-reference bases, so it is equally capable of representing 
/// substitutions, insertions, deletions, and complex indels. 
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub struct Variant {
    // leftmost position when n_tgt_bases > 0, or the position preceding an insertion
    pub tgt_pos0: SeqPos0,
    // for substitutions and deletions, the number of expected bases replaced by alt_bases         
    pub n_tgt_bases: ChromLength,
    // for substitutions and insertions, the bases replacing the expected bases    
    pub alt_bases: Option<UppercaseACGTN>,
    // whether tgt_pos0 is relative to a reference chromosome or a haplotype consensus
    #[serde(serialize_with = "serialize_haplotype")]
    pub haplotype: Haplotype,
}
impl Variant {
    /// Create a new Variant instance with the specified fields.
    pub fn new(
        tgt_pos0:    SeqPos0, 
        n_tgt_bases: ChromLength, 
        alt_bases:   &str,
        haplotype:   Haplotype,
    ) -> Self {
        Variant {
            tgt_pos0,
            n_tgt_bases,
            alt_bases: if alt_bases.is_empty() {
                None 
            } else { 
                Some(alt_bases.to_string()) 
            },
            haplotype,
        }
    }

    /// Get the signed difference in ref vs. alt length.
    pub fn alt_minus_ref(&self) -> i32 {
        self.alt_bases.as_ref().map_or(0, |alt| alt.len() as i32) - 
        self.n_tgt_bases as i32
    }

    /// Pack a variant into a string representation for printing to 
    /// variant_reads file.
    pub fn to_string(
        &self, 
        tgt_start0: SeqPos0, 
        avg_qual:   PhredQual
    ) -> String {
        let alt_bases = self.alt_bases
            .as_deref()
            .unwrap_or("-");
        format!(
            "{}:{}:{}:{}", 
            tgt_start0 + self.tgt_pos0,
            self.n_tgt_bases,
            alt_bases,
            avg_qual
        )
    }
}

//! SV/junction matching during benchmarking.

// imports
// use std::error::Error;
use std::cmp::Ordering;
use std::str::from_utf8_unchecked;
use serde::Serialize;
use rust_htslib::bcf;
use genomex::sam::junction::JunctionType;
use genomex::genome::{Chroms, TargetRegions};
use crate::junctions::FinalJunction;
use super::*;

// -----------------------------------------------------------------------------
// Handle mixtures of composite and non-composite genomes
// -----------------------------------------------------------------------------
/// Convert a HiFiRe3 chrom index into a non-composite canonical chrom name.
pub fn hf3_chrom_index_to_name(
    chrom_index1: u8, 
    is_composite_genome: bool,
    chroms: &Chroms, 
) -> Option<String> {
    simplify_chrom_name(&chroms.rev_index[&chrom_index1], is_composite_genome)
}
/// Convert a chrom name into a non-composite canonical chrom name.
pub fn simplify_chrom_name(
    chrom: &str, 
    is_composite_genome: bool
) -> Option<String> {
    let chrom = if is_composite_genome {
        let parts: Vec<_> = chrom.split("_").collect();
        let last_i0 = parts.len() - 1;
        if parts[last_i0] == "dm6" { return None } // drop fly composite alignments
        parts[0..last_i0].join("_") // drop the genome suffix
    } else {
        chrom.to_string()
    };
    if chrom.contains("_") { 
        None // drop non-canonical chroms 
    } else { 
        Some(chrom) 
    }
}

/// Check whether an SV breakpoint is on target.
pub fn sv_is_on_target(
    chrom_1: &str,
    pos1_1:  u32,
    chrom_2: &str,
    pos1_2:  u32,
    targets: &TargetRegions
) -> bool {
    !targets.has_data || 
    targets.get_pos_region(chrom_1, pos1_1).2 == INSIDE_REGION ||
    targets.get_pos_region(chrom_2, pos1_2).2 == INSIDE_REGION
}

// -----------------------------------------------------------------------------
// SV sources and types for HF3 vs. VCF comparison
// -----------------------------------------------------------------------------
/// Enum for tracking SV/junction sources.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Serialize)]
pub enum SvSource {
    #[serde(rename = "HF3")]
    HiFiRe3,
    #[serde(rename = "SNF")]
    Sniffles,
    #[serde(rename = "CSV")]
    CuteSV,
}

/// Enum for normalizing SV/junction types.
#[derive(PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SvType {
    Del,
    Dup,
    Inv,
    Trans,
    Bnd,
    Ins,
    Unknown,
}
impl SvType {
    pub fn from_vcf(jxn: &bcf::Record) ->  (Self, Vec<u8>) {
        let sv_type_vcf = get_vcf_sv_type_str(jxn);
        let sv_type = match sv_type_vcf.as_str() {
            "DEL" => SvType::Del,
            "DUP" => SvType::Dup,
            "INV" => SvType::Inv,
            "INS" => SvType::Ins,
            "BND" => match jxn.info(b"CHR2").string() {
                Ok(Some(chr2)) => unsafe { 
                    let (chrom_1, _) = get_vcf_left_node(&jxn)
                        .expect("Failed chrom in SvType::from_vcf");
                    let chrom_2 = from_utf8_unchecked(chr2[0]);
                    if chrom_1 == chrom_2 {
                        let alt  = get_vcf_alt(&jxn).expect("Failed alt in SvType::from_vcf");
                        match alt {
                            b"<DEL>" => SvType::Del,
                            b"<DUP>" => SvType::Dup,
                            b"<INV>" => SvType::Inv,
                            b"<INS>" => SvType::Ins,
                            _ => SvType::Bnd,
                        }
                    } else {
                        SvType::Trans
                    }
                },
                _ => SvType::Bnd,
            },
            _ => SvType::Unknown
        };
        (sv_type, sv_type_vcf.as_bytes().to_vec())
    }
    pub fn from_final_junction(item: &FinalJunction) -> Self {
        match item.jxn_type {
            x if x == JunctionType::Deletion as u8 => SvType::Del,
            x if x == JunctionType::Duplication as u8 => SvType::Dup,
            x if x == JunctionType::Inversion as u8 => SvType::Inv,
            x if x == JunctionType::Translocation as u8 => SvType::Trans,
            _ => SvType::Unknown
        }
    }
}

// -----------------------------------------------------------------------------
// SV junction calls
// -----------------------------------------------------------------------------
/// Table for aggregating SV/junction calls.
#[derive(Serialize)]
pub struct JunctionCall {
    pub sv_type: SvType,

    pub chrom_1: String,
    pub pos1_1:  u32,
    pub chrom_2: String,
    pub pos1_2:  u32,
    pub strands: String,

    pub source:  SvSource,

    pub sv_size:     u32,
    pub frag_size:   i32, // sign records whether it passed size filters
    pub stem_length: u32,

    pub offset:  i16,
    pub ins_seq: String,

    pub n_observed: u16,
    pub mapq: u8,
    pub divergence: f32,

    pub jxn_fail: u8,
    pub aln_fail: u8,
    pub all_chimeric: u8,

    pub mosaic: bool,
    pub vaf: f32,
    
    pub qnames: String,
}
impl JunctionCall {
    fn sort_key(
        &self,
    ) -> (&SvType, &str, u32, &str, u32, &str, &SvSource) {
        (
            &self.sv_type,
            self.chrom_1.as_str(),
            self.pos1_1,
            self.chrom_2.as_str(),
            self.pos1_2,
            self.strands.as_str(),
            &self.source,
        )
    }
}
impl PartialEq for JunctionCall {
    fn eq(&self, other: &Self) -> bool {
        self.sort_key() == other.sort_key()
    }
}
impl Eq for JunctionCall {}
impl PartialOrd for JunctionCall {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for JunctionCall {
    fn cmp(&self, other: &Self) -> Ordering {
        self.sort_key().cmp(&other.sort_key())
    }
}
// is_excluded_1:  u8,
// is_excluded_2:  u8,

// -----------------------------------------------------------------------------
//sorted QNAME lists
// -----------------------------------------------------------------------------
/// Assemble a sorted qnames concatenation from a comma-delimited String (HF3).
pub fn get_sorted_qnames_str(qnames: &str) -> String {
    let mut qnames: Vec<_> = qnames
        .split(",")
        .map(str::trim)
        .filter(|&s| !s.is_empty())
        .collect();
    sort_and_concatenate_qnames(&mut qnames[1..]) // since formatted as `,xxx,...`
}
/// Assemble a sorted qnames concatenation from a Vec<Vec<u8>> (VCF).
pub fn get_sorted_qnames_bytes(qnames: Vec<Vec<u8>>) -> String {
    let mut qnames: Vec<_> = qnames.iter().map(|x|{
        unsafe{ from_utf8_unchecked(x) }
    }).collect();
    sort_and_concatenate_qnames(&mut qnames)
}
fn sort_and_concatenate_qnames(qnames: &mut [&str]) -> String {
    qnames.sort_unstable();
    qnames.join(",")
}

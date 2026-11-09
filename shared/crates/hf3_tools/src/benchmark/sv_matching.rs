//! SV/junction matching during benchmarking.

// imports
// use std::error::Error;
use std::str::from_utf8_unchecked;
use serde::Serialize;
use rust_htslib::bcf;
use genomex::sam::junction::JunctionType;
use genomex::genome::Chroms;
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

// -----------------------------------------------------------------------------
// SV sources and types for HF3 vs. VCF comparison
// -----------------------------------------------------------------------------
/// Enum for tracking SV/junction sources.
#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SvSource {
    HiFiRe3,
    TruthSet,
}

/// Enum for normalizing SV/junction types.
#[derive(Serialize)]
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
    pub fn from_vcf(sv: &bcf::Record) -> Self {
        match sv.info(b"SVTYPE").string() {
            Ok(Some(tag)) => match tag[0] {
                b"DEL" => SvType::Del,
                b"DUP" => SvType::Dup,
                b"INV" => SvType::Inv,
                b"INS" => SvType::Ins,
                b"BND" => match sv.info(b"CHR2").string() {
                    Ok(Some(chr2)) => unsafe { 
                        let (chrom_1, _) = get_vcf_left_node(&sv)
                            .expect("Failed chrom in SvType::from_vcf");
                        let chrom_2 = from_utf8_unchecked(chr2[0]);
                        if chrom_1 == chrom_2 {
                            let alt  = get_vcf_alt(&sv).expect("Failed alt in SvType::from_vcf");
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
            },
            _ => SvType::Unknown
        }
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
    pub source:  SvSource,
    pub sv_type: SvType,
    pub sv_size: u32,
    pub n_observed: u16,
    pub chrom_1: String,
    pub pos1_1:  u32,
    pub chrom_2: String,
    pub pos1_2:  u32,
    pub jxn_seq: String,
    pub frag_size: i32,
    pub qnames: String,
}
    
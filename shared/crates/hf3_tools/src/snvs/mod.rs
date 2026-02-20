//! Module to handle SNV/indel analysis according to expectations
//! and encoding of HiFiRe3 libraries. 

// modules 
pub mod tags;
pub mod pileup;
pub mod variant;

// re-exports
pub use pileup::*;
pub use variant::*;

// dependencies
use std::error::Error;
use mdi::pub_key_constants;
use mdi::workflow::Config;
use genomex::genome::{Chroms, TargetRegions, Exclusions};

// constants
pub_key_constants!(
    // from environment variables
    SEQUENCING_PLATFORM
    LIBRARY_TYPE
);

/// Ensure that PacBio SNV analysis is performed on a library from 
/// the PacBioStrand sequencing platform.
pub fn check_pacbio_strand(tool: &str, cfg: &mut Config) -> Result<(), Box<dyn Error>> {
    cfg.set_string_env(&[SEQUENCING_PLATFORM, LIBRARY_TYPE]);
    let sequencing_platform = cfg.get_string(SEQUENCING_PLATFORM);
    let library_type        = cfg.get_string(LIBRARY_TYPE);
    if sequencing_platform != "PacBioStrand" || 
       library_type        != "HiFi" {
        return Err(format!(
            "{} requires PacBioStrand HiFi reads; found {} {} reads", 
            tool, sequencing_platform, library_type).into()
        );
    }
    Ok(())
}

/// The SnvAnalysisTool collects structs and methods for SNV analysis
/// at the genome level after all chromosomes have been processed.
pub struct SnvAnalysisTool {

    // global configuration parameters
    pub n_cpu:     u32,

    // chromosomes and regions
    pub chroms:     Chroms,
    pub targets:    TargetRegions,
    pub exclusions: Exclusions,
}

/// The SnvChromWorker tool collects structs and methods 
/// for SNV analysis while processing a single chromosome.
pub struct SnvChromWorker{

    // chromosome parsing
    pub chrom:       String,
    pub chrom_index: u8,

    // flag for the level of reads to include
    pub include_all_reads: bool,
    pub min_n_passes:      u8,

    // data structures for chromosome processing
    pub pileup:   ChromPileup,
    pub variants: ChromVariants,

    // output file paths
    pub pileup_file_path:   String,
    pub variants_file_path: String,
}

// SnvChromWorkerData enum allows difference types of metadata to be
// trasmitted to the main thread for aggregation.
pub enum SnvChromWorkerData {
    TotalAlnCount(usize),
    ErrorCorrectedAlignmentCount((String, usize)),
    PileupMetadata(PileupMetadata),
    VariantMetadata(VariantMetadata),
}

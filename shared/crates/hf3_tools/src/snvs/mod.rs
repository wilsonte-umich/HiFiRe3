//! Module to handle SNV/indel analysis according to expectations and encoding 
//! of HiFiRe3 libraries. 

// modules 
pub mod fragment;
pub mod tags;
pub mod simple_repeat;
pub mod variant;
pub mod consensus;
pub mod encoding;

// re-exports
pub use fragment::*;
pub use simple_repeat::*;
pub use variant::*;
// pub use encoding::*;

// dependencies
use std::error::Error;
use faimm::IndexedFasta;
use mdi::pub_key_constants;
use mdi::workflow::Config;
use genomex::genome::{Chroms, TargetRegions, Exclusions};
use genomex::sequence::Aligner;

/// type aliases
type SampleBit = u32;
type ChromName = String;
type ChromIndex = u8;
type ChromPos0 = u32;
type ChromPos1 = u32; // as used for 0-based _exclusive_ end positions (same as 1-based ends)
type ReadIndex = usize;
type VariantBases = String;

// constants
pub_key_constants!(
    // from environment variables
    SEQUENCING_PLATFORM
    LIBRARY_TYPE
);

/// Ensure that PacBio SNV analysis is performed on a library from the 
/// PacBioStrand sequencing platform.
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

/// SnvAnalysisTool collects tools for SNV analysis shared with all chromosome 
/// workers.
pub struct SnvAnalysisTool {

    // global configuration parameters
    pub n_cpu:     u32,

    // chromosomes and regions
    pub chroms:     Chroms,
    pub targets:    TargetRegions,
    pub exclusions: Exclusions,
    pub fa:         IndexedFasta,
}

/// SnvChromWorker collects tools for SNV analysis that are specific to 
/// processing of a single specific chromosome in parallel.
pub struct SnvChromWorker{

    // chromosome parsing
    pub chrom:       ChromName,
    pub chrom_index: ChromIndex,
    pub chrom_tid:   usize,
    pub simple_repeats: SimpleRepeats,

    // data structures for chromosome processing
    pub aligner: Aligner,
    // pub encodings: ReadEncodings,

    // output file paths
    pub variants_file_path:  String,
    pub encodings_file_path: String,
}

// SnvChromWorkerData enum allows difference types of metadata to be trasmitted 
// to the main thread for aggregation.
pub enum SnvChromWorkerData {
    TotalAlnCount(usize),
    UsableAlnCount((ChromName, usize)),
    // VariantMetadata(VariantMetadata),
    // EncodingMetadata(EncodingMetadata),
}

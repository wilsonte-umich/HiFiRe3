//! Module to handle SNV/indel analysis according to expectations and encoding 
//! of HiFiRe3 libraries. 

// modules 
pub mod fragment;
pub mod tags;
pub mod simple_repeat;
pub mod variant;
pub mod analyze_reads;
pub mod haplotype;
pub mod encoding;

// re-exports
pub use fragment::*;
pub use simple_repeat::*;
pub use variant::*;
pub use haplotype::*;
pub use encoding::*;

// dependencies
use std::error::Error;
use faimm::IndexedFasta;
use serde::Serialize;
use mdi::pub_key_constants;
use mdi::workflow::Config;
use genomex::genome::{Chroms, TargetRegions, Exclusions};
use genomex::sequence::Aligner;
use crate::tools::type_aliases::*;

// constants
pub_key_constants!(
    // from environment variables
    SEQUENCING_PLATFORM
    LIBRARY_TYPE
);
pub const MAX_HETEROZYGOUS_ZYGOSITY: f64 = 0.8;

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
    pub n_cpu: u32,

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
    pub chrom_index: ChromIndex1,
    pub chrom_tid:   usize,
    pub simple_repeats: SimpleRepeats,

    // processing tools
    pub aligner: Aligner,

    // metadata aggregation
    pub variant_tally:      VariantsTally,
    pub reads_on_reference: AlignmentEncodings,
    pub reads_on_haplotype: AlignmentEncodings,

    // output file paths
    pub variants_file_path:      String,
    pub reads_on_reference_path: String,
    pub reads_on_haplotype_path: String,
}

/// SnvChromWorkerData allows difference types of metadata to be trasmitted to 
/// the main thread for aggregation over the entire input.
pub enum SnvChromWorkerData {
    TotalAlnCount(usize),
    UsableAlnCount((ChromName, usize)),
    VariantMetadata(VariantMetadata),
    ClonalEncodingMetadata(EncodingMetadata),
    SubclonalEncodingMetadata(EncodingMetadata),
}

//! Count unique SNVs/indels in alignments and create a pileup.

// modules
mod chrom_worker;

// dependencies
use std::error::Error;
use crossbeam::channel::{bounded, unbounded};
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use genomex::genome::{Chroms, TargetRegions, Exclusions};
use crate::snvs::*;

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "analyze_snvs";
pub_key_constants!(
    // from environment variables
    N_CPU
    ANALYSIS_CHROMS_FILE
    // counter keys
    N_TOTAL_ALNS // split_by_chrom restricted input to on-target reads
    N_ALNS
    N_ALNS_BY_CHROM
    //-----------------------
    PILEUP_N_CHUNKS
    PILEUP_N_REPORTED_CHUNKS
    PILEUP_N_REPORTED_BASES
    PILEUP_REPORTED_COVERAGE
    PILEUP_N_VARIANT_BASES
    //-----------------------
    VARIANT_N_VARIANTS
    VARIANT_N_SUBSTITUTIONS
    VARIANT_N_INSERTIONS
    VARIANT_N_DELETIONS
    VARIANT_COVERAGE
);
const CHANNEL_CAPACITY: usize = 100;

// function called by hf3_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_usize_env( &[N_CPU]);
    cfg.set_string_env(&[ANALYSIS_CHROMS_FILE]);
                              
    // validate we are working with the expected read data type
    check_pacbio_strand(TOOL, &mut cfg)?;

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_TOTAL_ALNS, "usable on-target alignments processed, including non-error-corrected"),
        (N_ALNS,       "usable on-target error-corrected alignments processed"),
        (PILEUP_N_CHUNKS,          "number of error-corrected pileup chunks"),
        (PILEUP_N_REPORTED_CHUNKS, "number of pileup chunks reported in output"),
        (PILEUP_N_REPORTED_BASES,  "number of bases in reported pileup chunks"),
        (PILEUP_REPORTED_COVERAGE, "total base coverage in reported pileup chunks"),
        (PILEUP_N_VARIANT_BASES,   "number of variant bases in reported pileup chunks"),
        (VARIANT_N_VARIANTS,       "number of error-corrected SNV/indel variants reported"),
        (VARIANT_N_SUBSTITUTIONS,  "number of equal-length substitution variants"),
        (VARIANT_N_INSERTIONS,     "number of insertion variants"),
        (VARIANT_N_DELETIONS,      "number of deletion variants"),
        (VARIANT_COVERAGE,         "total base coverage at SNV/indel variants"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_ALNS_BY_CHROM,    "number of error-corrected alignments by on-target chromosome"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();

    // collect the working chromosomes
    let chroms = Chroms::new(&mut w.cfg);
    let targets = TargetRegions::from_env(&mut w, false);
    let on_target_chroms = targets.get_region_chroms(&chroms);
    chroms.write_chroms_file(w.cfg.get_string(ANALYSIS_CHROMS_FILE))?;

    // create the SNV analysis tool
    let tool = SnvAnalysisTool {
        n_cpu:      *w.cfg.get_usize(N_CPU) as u32,
        chroms:     chroms,
        targets:    targets,
        exclusions: Exclusions::from_env(&mut w, false),
    };

    // create channels for parallel processing
    let (tx_chrom, rx_chrom)    
        = unbounded::<(String, u8)>();
    let (tx_data, rx_data) 
        = bounded::<SnvChromWorkerData>(CHANNEL_CAPACITY);

    // spawn chromosome worker threads
    w.log.print("analyzing reads by chromosome");
    crossbeam::scope(|scope| {

        // workers: process one chromosome at a time
        let n_worker_threads = *w.cfg.get_usize(N_CPU) - 1; // leave one thread for collectors
        for _ in 0..n_worker_threads.max(1) {
            let rx_chrom = rx_chrom.clone();
            let tx_data = tx_data.clone();
            scope.spawn(|_| {
                chrom_worker::process_chrom(
                    &tool,
                    rx_chrom,
                    tx_data,
                ).unwrap();
            });
        }
        drop(tx_data);

        // transmit the chromosomes to be processed
        for (chrom, chrom_index) in on_target_chroms {
            tx_chrom.send((chrom, chrom_index)).unwrap();
        }
        drop(tx_chrom); 

        // collect metadata from chrom workers
        for metadata in rx_data {
            match metadata {
                SnvChromWorkerData::TotalAlnCount(count) => {
                    w.ctrs.add_to(N_TOTAL_ALNS, count);
                },
                SnvChromWorkerData::ErrorCorrectedAlignmentCount((chrom_name, count)) => {
                    w.ctrs.add_to(N_ALNS, count);
                    w.ctrs.add_to_keyed(N_ALNS_BY_CHROM, &chrom_name, count);
                },
                SnvChromWorkerData::PileupMetadata(md) => {
                    w.ctrs.add_to(PILEUP_N_CHUNKS,          md.n_chunks);
                    w.ctrs.add_to(PILEUP_N_REPORTED_CHUNKS, md.n_reported_chunks);
                    w.ctrs.add_to(PILEUP_N_REPORTED_BASES,  md.n_reported_bases);
                    w.ctrs.add_to(PILEUP_REPORTED_COVERAGE, md.reported_coverage);
                    w.ctrs.add_to(PILEUP_N_VARIANT_BASES,   md.n_variant_bases);
                },
                SnvChromWorkerData::VariantMetadata(md) => {
                    w.ctrs.add_to(VARIANT_N_VARIANTS,      md.n_variants);
                    w.ctrs.add_to(VARIANT_N_SUBSTITUTIONS, md.n_substitutions);
                    w.ctrs.add_to(VARIANT_N_INSERTIONS,    md.n_insertions);
                    w.ctrs.add_to(VARIANT_N_DELETIONS,     md.n_deletions);
                    w.ctrs.add_to(VARIANT_COVERAGE,        md.variant_coverage);
                },
            }
        }
    }).expect("Crossbeam scope panicked");

    // print counts
    w.ctrs.print_grouped(&[
        &[N_TOTAL_ALNS, N_ALNS],
        &[N_ALNS_BY_CHROM],
        &[PILEUP_N_CHUNKS, PILEUP_N_REPORTED_CHUNKS, PILEUP_N_REPORTED_BASES, PILEUP_REPORTED_COVERAGE, PILEUP_N_VARIANT_BASES],
        &[VARIANT_N_VARIANTS, VARIANT_N_SUBSTITUTIONS, VARIANT_N_INSERTIONS, VARIANT_N_DELETIONS, VARIANT_COVERAGE],
    ]);
    Ok(())
}

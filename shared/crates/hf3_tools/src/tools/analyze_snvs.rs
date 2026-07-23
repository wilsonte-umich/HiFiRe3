//! Count unique SNVs/indels in alignments and create a pileup.

// modules
mod chrom_worker;

// dependencies
use std::error::Error;
use crossbeam::channel::{bounded, unbounded};
use faimm::IndexedFasta;
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
    GENOME_FASTA
    // counter keys
    N_TOTAL_ALNS // split_by_chrom restricted input to on-target reads
    N_ALNS
    N_ALNS_BY_CHROM
    //-----------------------
    VARIANT_N_VARIANTS
    VARIANT_N_SUBSTITUTIONS
    VARIANT_N_INSERTIONS
    VARIANT_N_DELETIONS
    VARIANT_COUNT
    VARIANT_COVERAGE
    VARIANT_IS_SIMPLE_REPEAT
    //-----------------------
    ENCODING_N_READS
    ENCODING_N_REVERSE
    ENCODING_N_SPANS
    ENCODING_N_REF_BASES
    ENCODING_N_READ_BASES
    ENCODING_N_MATCH
    ENCODING_N_ALT
    ENCODING_N_MASKED
    ENCODING_N_DEL
    ENCODING_N_INS
    ENCODING_N_ALT_HIGH_QUAL
);
const CHANNEL_CAPACITY: usize = 100;

// function called by hf3_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_usize_env( &[N_CPU]);
    cfg.set_string_env(&[ANALYSIS_CHROMS_FILE, GENOME_FASTA]);
                              
    // validate we are working with the expected read data type
    check_pacbio_strand(TOOL, &mut cfg)?;

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_TOTAL_ALNS, "usable on-target alignments processed, including non-error-corrected"),
        (N_ALNS,       "usable on-target error-corrected alignments processed"),

        (VARIANT_N_VARIANTS,       "number of unique error-corrected SNV/indel variant types reported"),
        (VARIANT_N_SUBSTITUTIONS,  "number of equal-length substitution variant types"),
        (VARIANT_N_INSERTIONS,     "number of insertion variant types"),
        (VARIANT_N_DELETIONS,      "number of deletion variant types"),
        (VARIANT_COUNT,            "summed variant read count at all SNV/indel variant index positions"),
        (VARIANT_COVERAGE,         "summed read coverage at all SNV/indel variant index positions"),
        (VARIANT_IS_SIMPLE_REPEAT, "number of variants rejected as overlapping a simple repeat >= 6 units"),

        (ENCODING_N_READS,         "number of error-corrected reads subjected to encoding"),
        (ENCODING_N_REVERSE,       "number of those reads that aligned to the bottom reference strand"),
        (ENCODING_N_SPANS,         "number of unique genome alignment spans found in encoded reads"),
        (ENCODING_N_REF_BASES,     "number of genome bases covered by encoding alignment spans"),
        (ENCODING_N_READ_BASES,    "number of reference bases in encoded reads (M and D operations)"),
        (ENCODING_N_MATCH,         "number of reference-matched bases in encoded reads"),
        (ENCODING_N_ALT,           "number of unmasked alternative bases in encoded reads"),
        (ENCODING_N_MASKED,        "number of masked bases in encoded reads"),
        (ENCODING_N_DEL,           "number of deleted bases in encoded reads"),
        (ENCODING_N_INS,           "number of insertion operations in encoded reads (not the base count)"),
        (ENCODING_N_ALT_HIGH_QUAL, "number of high-quality alternative bases in encoded reads"),
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
    let genome_fasta = w.cfg.get_string(GENOME_FASTA).to_string();
    let tool = SnvAnalysisTool {
        n_cpu:      *w.cfg.get_usize(N_CPU) as u32,
        chroms:     chroms,
        targets:    targets,
        exclusions: Exclusions::from_env(&mut w, false),
        fa: IndexedFasta::from_file(&genome_fasta).expect("Error opening genome FASTA file")
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
                SnvChromWorkerData::UsableAlnCount((chrom_name, count)) => {
                    w.ctrs.add_to(N_ALNS, count);
                    w.ctrs.add_to_keyed(N_ALNS_BY_CHROM, &chrom_name, count);
                },
                // SnvChromWorkerData::VariantMetadata(md) => {
                //     w.ctrs.add_to(VARIANT_N_VARIANTS,      md.n_variants);
                //     w.ctrs.add_to(VARIANT_N_SUBSTITUTIONS, md.n_substitutions);
                //     w.ctrs.add_to(VARIANT_N_INSERTIONS,    md.n_insertions);
                //     w.ctrs.add_to(VARIANT_N_DELETIONS,     md.n_deletions);
                //     w.ctrs.add_to(VARIANT_COUNT,           md.variant_count);
                //     w.ctrs.add_to(VARIANT_COVERAGE,        md.variant_coverage);
                //     w.ctrs.add_to(VARIANT_IS_SIMPLE_REPEAT,md.is_simple_repeat);
                // },
                // SnvChromWorkerData::EncodingMetadata(md) => {
                //     w.ctrs.add_to(ENCODING_N_READS,        md.n_reads);
                //     w.ctrs.add_to(ENCODING_N_REVERSE,      md.n_reverse);
                //     w.ctrs.add_to(ENCODING_N_SPANS,        md.n_unique_spans);
                //     w.ctrs.add_to(ENCODING_N_REF_BASES,    md.n_ref_bases);
                //     w.ctrs.add_to(ENCODING_N_READ_BASES,   md.n_read_bases);
                //     w.ctrs.add_to(ENCODING_N_MATCH,        md.n_match);
                //     w.ctrs.add_to(ENCODING_N_ALT,          md.n_alt);
                //     w.ctrs.add_to(ENCODING_N_MASKED,       md.n_masked);
                //     w.ctrs.add_to(ENCODING_N_DEL,          md.n_del);
                //     w.ctrs.add_to(ENCODING_N_INS,          md.n_ins);
                //     w.ctrs.add_to(ENCODING_N_ALT_HIGH_QUAL,md.n_alt_high_qual);
                // },
            }
        }
    }).expect("Crossbeam scope panicked");

    // print counts
    w.ctrs.print_grouped(&[
        &[N_TOTAL_ALNS, N_ALNS],
        &[N_ALNS_BY_CHROM],
        &[
            VARIANT_N_VARIANTS, 
            VARIANT_N_SUBSTITUTIONS, 
            VARIANT_N_INSERTIONS, 
            VARIANT_N_DELETIONS, 
            VARIANT_COUNT,
            VARIANT_COVERAGE,
            VARIANT_IS_SIMPLE_REPEAT,
        ],
        &[
            ENCODING_N_READS,
            ENCODING_N_REVERSE,
            ENCODING_N_SPANS,
            ENCODING_N_REF_BASES,
            ENCODING_N_READ_BASES,
            ENCODING_N_MATCH,
            ENCODING_N_ALT,
            ENCODING_N_MASKED,
            ENCODING_N_DEL,
            ENCODING_N_INS,
            ENCODING_N_ALT_HIGH_QUAL,
        ],
    ]);
    Ok(())
}

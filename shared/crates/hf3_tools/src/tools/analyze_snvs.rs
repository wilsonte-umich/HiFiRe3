//! Analyze PacBio RE fragment haplotypes for clonal and subclonal variants.

// modules
mod chrom_worker;

// imports
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
    VARIANT_N_VARIANTS // variant tally
    VARIANT_N_SUBSTITUTIONS
    VARIANT_N_INSERTIONS
    VARIANT_N_DELETIONS
    VARIANT_N_CLONAL
    VARIANT_N_SUBCLONAL
    VARIANT_COUNT
    VARIANT_COVERAGE
    //-----------------------
    VARIANT_READS_N_READS // tally of reads with subclonal variants
    VARIANT_READS_N_INDEL_ONLY
    VARIANT_READS_N_ONE_SNV
    VARIANT_READS_N_TWO_SNV
    VARIANT_READS_N_THREE_SNV
    VARIANT_READS_N_FOUR_SNV
    VARIANT_READS_N_FIVE_SNV
    //-----------------------
    ON_REFERENCE_SPANS // read encodings
    ON_REFERENCE_READS
    ON_REFERENCE_REF_BASES
    ON_REFERENCE_READ_BASES
    ON_REFERENCE_MATCH
    ON_REFERENCE_ALT
    ON_REFERENCE_DEL
    ON_REFERENCE_INS
    ON_REFERENCE_MASKED
    //-----------------------
    ON_HAPLOTYPE_SPANS
    ON_HAPLOTYPE_READS
    ON_HAPLOTYPE_REF_BASES
    ON_HAPLOTYPE_READ_BASES
    ON_HAPLOTYPE_MATCH
    ON_HAPLOTYPE_ALT
    ON_HAPLOTYPE_DEL
    ON_HAPLOTYPE_INS
    ON_HAPLOTYPE_MASKED
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

        (VARIANT_N_VARIANTS,        "number of unique SNV/indel variants reported"),
        (VARIANT_N_SUBSTITUTIONS,   "number of equal-length substitution variants"),
        (VARIANT_N_INSERTIONS,      "number of insertion variants"),
        (VARIANT_N_DELETIONS,       "number of deletion variants"),
        (VARIANT_N_CLONAL,          "number of clonal variants"),
        (VARIANT_N_SUBCLONAL,       "number of subclonal variants"),
        (VARIANT_COUNT,             "summed variant read count at all index positions"),
        (VARIANT_COVERAGE,          "summed read coverage at all SNV/indel index positions"),

        (VARIANT_READS_N_READS,     "number of reads with at least one subclonal variant reported"),
        (VARIANT_READS_N_INDEL_ONLY,"number of reads with only subclonal indels reported"),
        (VARIANT_READS_N_ONE_SNV,   "number of reads with one subclonal SNV reported"),
        (VARIANT_READS_N_TWO_SNV,   "number of reads with two subclonal SNVs reported"),
        (VARIANT_READS_N_THREE_SNV, "number of reads with three subclonal SNVs reported"),
        (VARIANT_READS_N_FOUR_SNV,  "number of reads with four subclonal SNVs reported"),
        (VARIANT_READS_N_FIVE_SNV,  "number of reads with five or more subclonal SNVs reported"),

        (ON_REFERENCE_SPANS,        "number of unique genome alignment spans found in clonal read_on_ref analysis"),
        (ON_REFERENCE_READS,        "number of error-corrected reads subjected to read_on_ref analysis"),
        (ON_REFERENCE_REF_BASES,    "number of genome bases covered by ON_REFERENCE_SPANS"),
        (ON_REFERENCE_READ_BASES,   "number of reference bases in ON_REFERENCE_READS (M and D operations)"),
        (ON_REFERENCE_MATCH,        "number of reference-matched bases in ON_REFERENCE_READS"),
        (ON_REFERENCE_ALT,          "number of unmasked alternative bases in ON_REFERENCE_READS"),
        (ON_REFERENCE_DEL,          "number of deleted bases in ON_REFERENCE_READS"),
        (ON_REFERENCE_INS,          "number of insertion operations in ON_REFERENCE_READS (not the base count)"),
        (ON_REFERENCE_MASKED,       "number of masked bases in ON_REFERENCE_READS"),

        (ON_HAPLOTYPE_SPANS,        "number of unique genome alignment spans found in subclonal read_on_hap analysis"),
        (ON_HAPLOTYPE_READS,        "number of error-corrected reads subjected to read_on_hap analysis"),
        (ON_HAPLOTYPE_REF_BASES,    "number of genome bases covered by ON_HAPLOTYPE_SPANS"),
        (ON_HAPLOTYPE_READ_BASES,   "number of reference bases in ON_HAPLOTYPE_READS (M and D operations)"),
        (ON_HAPLOTYPE_MATCH,        "number of reference-matched bases in ON_HAPLOTYPE_READS"),
        (ON_HAPLOTYPE_ALT,          "number of unmasked alternative bases in ON_HAPLOTYPE_READS"),
        (ON_HAPLOTYPE_DEL,          "number of deleted bases in ON_HAPLOTYPE_READS"),
        (ON_HAPLOTYPE_INS,          "number of insertion operations in ON_HAPLOTYPE_READS (not the base count)"),
        (ON_HAPLOTYPE_MASKED,       "number of masked bases in ON_HAPLOTYPE_READS"),
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
        let n_worker_threads = *w.cfg.get_usize(N_CPU);
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
                SnvChromWorkerData::VariantMetadata(md) => {
                    w.ctrs.add_to(VARIANT_N_VARIANTS,      md.n_variants);
                    w.ctrs.add_to(VARIANT_N_SUBSTITUTIONS, md.n_substitutions);
                    w.ctrs.add_to(VARIANT_N_INSERTIONS,    md.n_insertions);
                    w.ctrs.add_to(VARIANT_N_DELETIONS,     md.n_deletions);
                    w.ctrs.add_to(VARIANT_N_CLONAL,        md.n_clonal);
                    w.ctrs.add_to(VARIANT_N_SUBCLONAL,     md.n_subclonal);
                    w.ctrs.add_to(VARIANT_COUNT,           md.variant_count);
                    w.ctrs.add_to(VARIANT_COVERAGE,        md.variant_coverage);
                },
                SnvChromWorkerData::VariantReadsMetadata(md) => {
                    w.ctrs.add_to(VARIANT_READS_N_READS,   md.n_variant_reads);
                    w.ctrs.add_to(VARIANT_READS_N_INDEL_ONLY, md.n_indel_only);
                    w.ctrs.add_to(VARIANT_READS_N_ONE_SNV, md.n_one_snv);
                    w.ctrs.add_to(VARIANT_READS_N_TWO_SNV, md.n_two_snv);
                    w.ctrs.add_to(VARIANT_READS_N_THREE_SNV, md.n_three_snv);
                    w.ctrs.add_to(VARIANT_READS_N_FOUR_SNV,md.n_four_snv);
                    w.ctrs.add_to(VARIANT_READS_N_FIVE_SNV,md.n_five_snv);
                },
                SnvChromWorkerData::ReadsOnReferenceMetadata(md) => {
                    w.ctrs.add_to(ON_REFERENCE_SPANS,     md.n_unique_spans);
                    w.ctrs.add_to(ON_REFERENCE_READS,     md.n_reads);
                    w.ctrs.add_to(ON_REFERENCE_REF_BASES, md.n_ref_bases);
                    w.ctrs.add_to(ON_REFERENCE_READ_BASES,md.n_read_bases);
                    w.ctrs.add_to(ON_REFERENCE_MATCH,     md.n_match);
                    w.ctrs.add_to(ON_REFERENCE_ALT,       md.n_alt);
                    w.ctrs.add_to(ON_REFERENCE_DEL,       md.n_del);
                    w.ctrs.add_to(ON_REFERENCE_INS,       md.n_ins);
                    w.ctrs.add_to(ON_REFERENCE_MASKED,    md.n_masked);
                },
                SnvChromWorkerData::ReadsOnHaplotypeMetadata(md) => {
                    w.ctrs.add_to(ON_HAPLOTYPE_SPANS,     md.n_unique_spans);
                    w.ctrs.add_to(ON_HAPLOTYPE_READS,     md.n_reads);
                    w.ctrs.add_to(ON_HAPLOTYPE_REF_BASES, md.n_ref_bases);
                    w.ctrs.add_to(ON_HAPLOTYPE_READ_BASES,md.n_read_bases);
                    w.ctrs.add_to(ON_HAPLOTYPE_MATCH,     md.n_match);
                    w.ctrs.add_to(ON_HAPLOTYPE_ALT,       md.n_alt);
                    w.ctrs.add_to(ON_HAPLOTYPE_DEL,       md.n_del);
                    w.ctrs.add_to(ON_HAPLOTYPE_INS,       md.n_ins);
                    w.ctrs.add_to(ON_HAPLOTYPE_MASKED,    md.n_masked);
                },
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
            VARIANT_N_CLONAL,
            VARIANT_N_SUBCLONAL,
            VARIANT_COUNT,
            VARIANT_COVERAGE,
        ],
        &[
            VARIANT_READS_N_READS,
            VARIANT_READS_N_INDEL_ONLY,
            VARIANT_READS_N_ONE_SNV,
            VARIANT_READS_N_TWO_SNV,
            VARIANT_READS_N_THREE_SNV,
            VARIANT_READS_N_FOUR_SNV,
            VARIANT_READS_N_FIVE_SNV,
        ],
        &[
            ON_REFERENCE_SPANS,
            ON_REFERENCE_READS,
            ON_REFERENCE_REF_BASES,
            ON_REFERENCE_READ_BASES,
            ON_REFERENCE_MATCH,
            ON_REFERENCE_ALT,
            ON_REFERENCE_DEL,
            ON_REFERENCE_INS,
            ON_REFERENCE_MASKED,
        ],
        &[
            ON_HAPLOTYPE_SPANS,
            ON_HAPLOTYPE_READS,
            ON_HAPLOTYPE_REF_BASES,
            ON_HAPLOTYPE_READ_BASES,
            ON_HAPLOTYPE_MATCH,
            ON_HAPLOTYPE_ALT,
            ON_HAPLOTYPE_DEL,
            ON_HAPLOTYPE_INS,
            ON_HAPLOTYPE_MASKED,
        ],
    ]);
    Ok(())
}

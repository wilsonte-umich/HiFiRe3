//! Extract read paths, junctions, and alignments from chrom-level 
//! BAM files and group SVs in preparation for structural variant 
//! analysis and visualization in the R Shiny app.

// modules
mod chrom_worker;

// dependencies
use std::error::Error;
use rustc_hash::FxHashMap;
use crossbeam::channel::{bounded, unbounded};
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use genomex::genome::{Chroms, TargetRegions, Genes, Exclusions};
use genomex::sam::junction::JunctionType;
use crate::junctions::*;

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "analyze_svs";
pub_key_constants!(
    // from environment variables
    N_CPU
    GROUP_BREAKPOINT_DISTANCE
    GROUP_STEM_DISTANCE
    DEDUPLICATE_READS
    DATA_NAME
    SEQUENCING_PLATFORM
    INDEX_FILE_PREFIX_WRK
    SV_READ_PATHS_FILE
    SV_ALIGNMENTS_FILE
    SV_COVERAGE_FILE
    SV_FINAL_JUNCTIONS_FILE_1
    SV_FINAL_JUNCTIONS_FILE_2
    ANALYSIS_CHROMS_FILE
    // counter keys
    N_READS // split_by_chrom restricts input to on-target reads
    N_SV_READS
    N_FINAL_JUNCTIONS
    N_FIRST_ALNS
    N_UNIQ_FIRST_ALNS
    N_DISTAL_ALNS
    N_UNIQ_DISTAL_ALNS
    N_READS_BY_CHROM
    N_JXNS_BY_OFFSET
    N_JXNS_BY_TYPE
    N_JXNS_BY_CHIMERICITY
    N_JXN_BY_GENOMICITY
    N_JXNS_BY_COUNT
    N_SINGLETONS_BY_CHIMERICITY
    N_VALIDATED_BY_CHIMERICITY
);
const CHANNEL_CAPACITY: usize = 100;
const JXN_CAPACITY:     usize = 100_000;

// function called by hf3_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_usize_env( &[N_CPU, GROUP_BREAKPOINT_DISTANCE, GROUP_STEM_DISTANCE]);
    cfg.set_bool_env(  &[DEDUPLICATE_READS]);
    cfg.set_string_env(&[DATA_NAME, SEQUENCING_PLATFORM, SV_READ_PATHS_FILE,
                              INDEX_FILE_PREFIX_WRK, SV_ALIGNMENTS_FILE, SV_COVERAGE_FILE,
                              SV_FINAL_JUNCTIONS_FILE_1, SV_FINAL_JUNCTIONS_FILE_2,
                              ANALYSIS_CHROMS_FILE]);

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_READS, "usable on-target reads processed"),
        (N_SV_READS,        "reads with committed SV junctions"),
        (N_FINAL_JUNCTIONS, "committed SV junctions"),
        (N_FIRST_ALNS,      "first alignments in reads"),
        (N_UNIQ_FIRST_ALNS, "unique first alignments in reads"),
        (N_DISTAL_ALNS,     "distal alignments in reads with SVs"),
        (N_UNIQ_DISTAL_ALNS,"unique distal alignments in reads with SVs"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_READS_BY_CHROM,    "number of reads by on-target chromosome"),
        (N_JXNS_BY_OFFSET,    "number of final junctions by breakpoint offset group"),
        (N_JXNS_BY_TYPE,      "number of final junctions by type"),
        (N_JXNS_BY_CHIMERICITY, "number of final junctions by chimericity"),
        (N_JXN_BY_GENOMICITY, "number of final junctions by genomicity"),
        (N_JXNS_BY_COUNT,     "number of final junctions by supporting read count"),
        (N_SINGLETONS_BY_CHIMERICITY, "number of singleton junctions (n==1)  by chimericity"),
        (N_VALIDATED_BY_CHIMERICITY,  "number of validated junctions (n >=3 )by chimericity"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();

    // collect the working chromosomes
    let chroms = Chroms::new(&mut w.cfg);
    let targets = TargetRegions::from_env(&mut w, false);
    let on_target_chroms = targets.get_region_chroms(&chroms);
    chroms.write_chroms_file(w.cfg.get_string(ANALYSIS_CHROMS_FILE))?;

    // create the junction analysis tool
    let mut tool = JunctionAnalysisTool {
        data_name:  w.cfg.get_string(DATA_NAME).to_string(),
        chroms:     chroms,
        targets:    targets,
        genes:      Genes::from_env(     &mut w, false, 0),
        exclusions: Exclusions::from_env(&mut w, false),
        group_breakpoint_distance: *w.cfg.get_usize(GROUP_BREAKPOINT_DISTANCE),
        group_stem_distance:       *w.cfg.get_usize(GROUP_STEM_DISTANCE) as u32,
        is_ont:                     w.cfg.get_string(SEQUENCING_PLATFORM).to_lowercase() == "ONT",
        deduplicate_reads:         *w.cfg.get_bool(DEDUPLICATE_READS),
        final_jxns_file_1:          w.cfg.get_string(SV_FINAL_JUNCTIONS_FILE_1).to_string(),
        final_jxns_file_2:          w.cfg.get_string(SV_FINAL_JUNCTIONS_FILE_2).to_string(),
    };

    // create channels for parallel processing
    let (tx_chrom, rx_chrom)    
        = unbounded::<(String, u8)>();
    let (tx_data, rx_data) 
        = bounded::<ChromWorkerData>(CHANNEL_CAPACITY);

    // create collectors to receive results from workers
    let mut read_paths:  Vec<SvReadPath>       = Vec::with_capacity(JXN_CAPACITY);
    let mut distal_alns: Vec<AlignmentSegment> = Vec::with_capacity(JXN_CAPACITY);
    let mut jxns: FxHashMap<OrderedJunction, JunctionInstances> = FxHashMap::default();

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
                    &tool.chroms,
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

        // collect read paths from chrom workers
        for data in rx_data {
            match data {
                ChromWorkerData::ReadPath(read_path) => 
                    read_paths.push(read_path),
                ChromWorkerData::DistalAln(distal_aln) => 
                    distal_alns.push(distal_aln),
                ChromWorkerData::Junction((ordered_jxn, jxn_instance)) => {
                    jxns.entry(ordered_jxn)
                        .or_insert_with(|| JunctionInstances::new())
                        .add_instance(jxn_instance);
                },
                ChromWorkerData::ChromReadCount((chrom_name, count)) => {
                    w.ctrs.add_to_keyed(N_READS_BY_CHROM, &chrom_name, count);
                    w.ctrs.add_to(N_READS, count);
                },
            }
        }
    }).expect("Crossbeam scope panicked");

    // sort and print aggregated SV read paths
    w.ctrs.add_to(N_SV_READS, read_paths.len());
    SvReadPath::write_sorted(
        read_paths, 
        w.cfg.get_string(SV_READ_PATHS_FILE)
    )?;

    // fuzzy group final junctions
    let mut final_jxns = fuzzy_match_junctions(jxns, &mut tool)?;
    w.ctrs.add_to(N_FINAL_JUNCTIONS, final_jxns.len());

    // merge distal alignments into first alignments
    // use AlignmentSegments to update final junction breakpoint coverage fields
    let (n_first_alns, n_uniq_first_alns, n_distal_alns, n_uniq_distal_alns)
         = AlignmentSegment::write_merged(
        distal_alns,
        w.cfg.get_string(INDEX_FILE_PREFIX_WRK),
             w.cfg.get_string(SV_ALIGNMENTS_FILE),
             w.cfg.get_string(SV_COVERAGE_FILE),
        &tool.chroms,
        &mut final_jxns,
    )?;
    w.ctrs.add_to(N_FIRST_ALNS,       n_first_alns);
    w.ctrs.add_to(N_UNIQ_FIRST_ALNS,  n_uniq_first_alns);
    w.ctrs.add_to(N_DISTAL_ALNS,      n_distal_alns);
    w.ctrs.add_to(N_UNIQ_DISTAL_ALNS, n_uniq_distal_alns);

    // collect final junction statistics
    for jxn in &final_jxns {
        let key = jxn.get_offset_type();
        w.ctrs.increment_keyed(N_JXNS_BY_OFFSET,      key);
        let key = JunctionType::str_from_u8(jxn.jxn_type);
        w.ctrs.increment_keyed(N_JXNS_BY_TYPE,        key);
        let key = jxn.get_chimericity();
        w.ctrs.increment_keyed(N_JXNS_BY_CHIMERICITY, key);
        if jxn.n_instances_dedup == 1 {
            w.ctrs.increment_keyed(N_SINGLETONS_BY_CHIMERICITY, key);
        } else if jxn.n_instances_dedup >= 3 {
            w.ctrs.increment_keyed(N_VALIDATED_BY_CHIMERICITY,  key);
        }        
        let key = jxn.get_genomicity();
        w.ctrs.increment_keyed(N_JXN_BY_GENOMICITY,   key);
        let key = jxn.get_coverage_type();
        w.ctrs.increment_keyed(N_JXNS_BY_COUNT,       key);
    }

    // write final junctions files in two sort orders
    FinalJunction::write_sorted(final_jxns, &tool);

    // print counts
    w.ctrs.print_grouped(&[
        &[N_READS, N_SV_READS],
        &[N_READS_BY_CHROM],
        &[N_FIRST_ALNS, N_UNIQ_FIRST_ALNS, N_DISTAL_ALNS, N_UNIQ_DISTAL_ALNS],
        &[N_FINAL_JUNCTIONS], 
        &[N_JXNS_BY_TYPE],
        &[N_JXNS_BY_OFFSET],
        &[N_JXNS_BY_COUNT],
        &[N_JXN_BY_GENOMICITY],
        &[N_JXNS_BY_CHIMERICITY],
        &[N_SINGLETONS_BY_CHIMERICITY],
        &[N_VALIDATED_BY_CHIMERICITY],
    ]);
    Ok(())
}

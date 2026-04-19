//! Merge SVs and supporting data from multiple input libraries of different types.

// dependencies
use std::error::Error;
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
    MERGE_INPUT_DIRS
    SV_READ_PATHS_FILE
    SV_ALIGNMENTS_FILE
    SV_COVERAGE_FILE // currently unfilled
    SV_FINAL_JUNCTIONS_FILE_1
    SV_FINAL_JUNCTIONS_FILE_2
    ANALYSIS_CHROMS_FILE
    SV_SAMPLES_FILE
    // counter keys
    N_FINAL_JUNCTIONS_IN
    N_FINAL_JUNCTIONS_OUT
    N_JXNS_BY_OFFSET
    N_JXNS_BY_TYPE
    N_JXNS_BY_CHIMERICITY
    N_JXN_BY_GENOMICITY
    N_JXNS_BY_COUNT
    N_SINGLETONS_BY_CHIMERICITY
    N_VALIDATED_BY_CHIMERICITY
);

// function called by hf3_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_usize_env( &[N_CPU]);
    cfg.set_string_env(&[SV_READ_PATHS_FILE, SV_ALIGNMENTS_FILE, SV_COVERAGE_FILE,
                              SV_FINAL_JUNCTIONS_FILE_1, SV_FINAL_JUNCTIONS_FILE_2,
                              ANALYSIS_CHROMS_FILE, SV_SAMPLES_FILE]);

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_FINAL_JUNCTIONS_IN, "input final junctions in all libraries"),
        (N_FINAL_JUNCTIONS_OUT,"output final junctions after merging libraries"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_JXNS_BY_OFFSET,      "number of final junctions by breakpoint offset group"),
        (N_JXNS_BY_TYPE,        "number of final junctions by type"),
        (N_JXNS_BY_CHIMERICITY, "number of final junctions by chimericity"),
        (N_JXN_BY_GENOMICITY,   "number of final junctions by genomicity"),
        (N_JXNS_BY_COUNT,       "number of final junctions by supporting read count"),
        (N_SINGLETONS_BY_CHIMERICITY, "number of singleton junctions (n==1)  by chimericity"),
        (N_VALIDATED_BY_CHIMERICITY,  "number of validated junctions (n >=3 )by chimericity"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();
    let n_cpu = *w.cfg.get_usize(N_CPU) as u32;

    // collect the working chromosomes
    let chroms = Chroms::new(&mut w.cfg);
    let targets = TargetRegions::from_env(&mut w, false);
    chroms.write_chroms_file(w.cfg.get_string(ANALYSIS_CHROMS_FILE))?;

    // parse the input directories
    let merge_input_dirs = w.cfg.get_string(MERGE_INPUT_DIRS).to_string();
    let merge_input_dirs = merge_input_dirs.split(',').collect::<Vec<&str>>();

    // create the junction analysis tool
    let tool = JunctionAnalysisTool {
        n_cpu:      *w.cfg.get_usize(N_CPU) as u32,
        chroms:     chroms,
        targets:    targets,
        genes:      Genes::from_env(     &mut w, false, 0),
        exclusions: Exclusions::from_env(&mut w, false),
        group_breakpoint_distance: 20, // the same values for most library types except ONT, which are usually more lenient
        group_stem_distance:       5,
        is_ont:                    false,
        deduplicate_reads:         false,
        final_jxns_file_1:         w.cfg.get_string(SV_FINAL_JUNCTIONS_FILE_1).to_string(),
        final_jxns_file_2:         w.cfg.get_string(SV_FINAL_JUNCTIONS_FILE_2).to_string(),
    };

    // merge samples files, bit shifting later samples as needed
    let jxns_in = Sample::parse_merge_samples(&mut w.cfg, &merge_input_dirs)?;
    let n_jxns_in = jxns_in.len();

    // sort and print merged SV read paths
    w.log.print("sorting and printing merged SV read paths");
    SvReadPath::merge_and_write_sorted(
        &merge_input_dirs, 
        w.cfg.get_string(SV_READ_PATHS_FILE),
        n_cpu,
    )?;

    // sort and print merged alignment segments
    w.log.print("sorting and printing merged alignment segments");
    AlignmentSegment::merge_and_write_sorted(
        &merge_input_dirs, 
        w.cfg.get_string(SV_ALIGNMENTS_FILE),
        n_cpu,
    )?;

    // fuzzy group final junctions
    w.log.print("fuzzy grouping final junctions");
    let final_jxns = fuzzy_merge_junctions(jxns_in, &tool)?;
    w.ctrs.add_to(N_FINAL_JUNCTIONS_IN,  n_jxns_in);
    w.ctrs.add_to(N_FINAL_JUNCTIONS_OUT, final_jxns.len());

    // collect final junction statistics
    w.log.print("collecting final junction statistics");
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
    w.log.print("sorting and writing final junctions");
    FinalJunction::write_sorted(final_jxns, &tool);

    // print counts
    w.ctrs.print_grouped(&[
        &[N_FINAL_JUNCTIONS_IN, N_FINAL_JUNCTIONS_OUT],
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

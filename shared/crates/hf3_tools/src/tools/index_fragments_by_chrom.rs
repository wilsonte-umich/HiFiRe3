//! Extract read paths, junctions, and alignments from chrom-level 
//! BAM files in preparation for structural variant analysis and 
//! visualization in the R Shiny app.

// modules
mod chrom_worker;

// dependencies
use std::error::Error;
use rustc_hash::FxHashMap;
use crossbeam::channel::{bounded, unbounded};
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use genomex::genome::{Chroms, TargetRegions};
use crate::junctions::*;

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "analyze_svs";
pub_key_constants!(
    // from environment variables
    N_CPU
    UNIQUE_JXNS_FILE
    SV_READ_PATHS_FILE
    DISTAL_ALNS_FILE
    // counter keys
    N_READS_ON_TARGET
    N_READS_BY_CHROM
);
const CHANNEL_CAPACITY: usize = 100;

// function called by hf3_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[UNIQUE_JXNS_FILE, SV_READ_PATHS_FILE, DISTAL_ALNS_FILE]);
    cfg.set_usize_env(&[N_CPU]);

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_READS_ON_TARGET, "on-target reads processed"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_READS_BY_CHROM, "number of reads by on-target chromosome"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.print("initializing");

    // collect the working chromosomes
    let chroms = Chroms::new(&mut w.cfg);
    let targets = TargetRegions::new(&mut w);
    let on_target_chroms = targets.get_on_target_chroms(&chroms);

    // create channels for fan-in/fan-out parallel processing
    let (tx_chrom, rx_chrom)    
        = unbounded::<(String, usize)>();
    let (tx_jxn, rx_jxn) 
        = bounded::<(OrderedJunction, JunctionInstance)>(CHANNEL_CAPACITY);
    let (tx_read_path, rx_read_path) 
        = bounded::<SvReadPath>(CHANNEL_CAPACITY);
    let (tx_distal_aln, rx_distal_aln) 
        = bounded::<AlignmentSegment>(CHANNEL_CAPACITY);

    // collectors: receive results from workers
    let mut read_paths: Vec<SvReadPath> = Vec::with_capacity(100000);
    let mut distal_alns: Vec<AlignmentSegment> = Vec::with_capacity(100000);
    let mut jxns: FxHashMap<OrderedJunction, JunctionInstances> = FxHashMap::default();

    // spawn threads
    crossbeam::scope(|scope| {

        // workers: process one chromosome at a time
        let n_worker_threads = *w.cfg.get_usize(N_CPU) - 1;
        for _ in 0..n_worker_threads {
            let rx_chrom  = rx_chrom.clone();
            let tx_jxn     = tx_jxn.clone();
            let tx_read_path = tx_read_path.clone();
            let tx_distal_aln = tx_distal_aln.clone();
            scope.spawn(|_| {
                chrom_worker::process_chrom(
                    &chroms,
                    rx_chrom,
                    tx_jxn,
                    tx_read_path,
                    tx_distal_aln,
                ).unwrap();
            });
        }

        // Collect read paths from chrom workers.
        for read_path in rx_read_path {
            read_paths.push(read_path);
        }

        // Collect distal alignments from chrom workers (first read alignment printed by workers).
        for distal_aln in rx_distal_aln {
            distal_alns.push(distal_aln);
        }

        // Collect junctions from chrom workers.
        for (ordered_jxn, jxn_instance) in rx_jxn {
            jxns.entry(ordered_jxn)
                .or_insert_with(|| JunctionInstances::new())
                .add_instance(jxn_instance);
        }

        // transmit the chromosomes to be processed
        for (chrom, chrom_index) in on_target_chroms {
            w.log.print(&format!("sending chrom {}", chrom));
            tx_chrom.send((chrom, chrom_index)).unwrap();
        }

        // close channels
        drop(tx_chrom); 
        drop(tx_jxn);
        drop(tx_read_path);
        drop(tx_distal_aln);

    }).expect("Crossbeam scope panicked");

    // sort and print read paths
    SvReadPath::write_sorted(
        read_paths, 
        w.cfg.get_string(SV_READ_PATHS_FILE)
    )?;

    // sort and print distal alignments (merged downstream with first alignments)
    AlignmentSegment::write_sorted_with_count(
        distal_alns, 
        w.cfg.get_string(DISTAL_ALNS_FILE)
    )?;

    // sort, fuzzy group and print junctions
    write_sorted_with_count(
        jxns, 
        w.cfg.get_string(UNIQUE_JXNS_FILE)
    )?;

    // print counts
    w.ctrs.print_grouped(&[
        &[N_READS_ON_TARGET],
        &[N_READS_BY_CHROM],
    ]);
    Ok(())
}

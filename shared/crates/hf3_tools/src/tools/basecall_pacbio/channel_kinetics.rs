//! Receiver to collect and print kinetics instances from worker threads.

// modules
pub mod kinetics;

// dependencies
use std::error::Error;
use std::fs::File;
use std::io::Write;
use rustc_hash::FxHashMap;
use flate2::write::GzEncoder;
use flate2::Compression;
use crossbeam::channel::Receiver;
use mdi::pub_key_constants;
use mdi::workflow::Config;
use super::KineticsInstance;
use kinetics::StrandMetadata;

// constants
pub_key_constants!(
    // from environment variables
    PACBIO_BASECALL_KINETICS // output kinetics file path
);
const MAX_CONTEXT_VALUES: usize = 1_000_000; // print up to one million entries per base context

/// Collect kinetics instances from workers. Write up to 1M frame count
/// trios for each base context type to PACBIO_BASECALL_KINETICS file.
pub fn collect_kinetics(
    rx_kinetics:  Receiver<KineticsInstance>,
) -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[PACBIO_BASECALL_KINETICS]);

    // initalize the kinetics output file
    let file = File::create(cfg.get_string(PACBIO_BASECALL_KINETICS))?;
    let mut file = GzEncoder::new(file, Compression::default());
    writeln!(
        file,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        "base_context",
        "is_heteroduplex",
        "ref_is_known",
        "is_ref",
        "ipd_before",
        "pulse_width",
        "ipd_after"
    )?;

    // remember how many instances have been printed by context type
    let mut context_counts: FxHashMap<StrandMetadata, usize> = FxHashMap::default();

    // process kinetics instances as they arrive
    for kinetics_instance in rx_kinetics.iter() {
        let md = &kinetics_instance.strand_metadata;
        let context_count = context_counts.get(md).unwrap_or(&0);
        if *context_count >= MAX_CONTEXT_VALUES {
            continue; // skip writing more entries for this context type
        }
        let fc = &kinetics_instance.frame_counts;
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            md.base_context.iter().map(|&b| b as char).collect::<String>(),
            md.is_heteroduplex as u8,
            md.ref_is_known    as u8,
            md.is_ref          as u8,
            fc.ip_before,
            fc.pw,
            fc.ip_after,
        )?;
        context_counts.insert(kinetics_instance.strand_metadata, context_count + 1);
    }
    Ok(())
}

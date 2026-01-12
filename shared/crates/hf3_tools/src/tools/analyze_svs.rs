//! This is a template for constructing a workflow tool using MDI framework components.
//! 
//! A workflow tool is called as `my_tools xxx` from the command line or a script and
//! performs a specific data processing task.
//! 
//! Replace this comment block with a description of the tool's purpose, actions
//! performed, expected inputs, and the generated outputs.

// modules
// mod sub_module; // declare sub-modules as needed, e.g., tool helpers, data formats, etc.

// dependencies
use std::error::Error;
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use mdi::RecordStreamer;
use genomex::sam::{SamRecord, flag};
use genomex::genome::Chroms;
use crate::formats::hf_tags::*;
use crate::junctions::JxnFailureFlag;

// Tool structure, for passing data workers to record processing functions.
struct Tool {
    my_worker: usize,
}

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "analyze_something"; // tool name for logging, typically in format verb_noun
pub_key_constants!(
    // from environment variables
    MY_OPTION 
    MY_FLAG
    // derived configuration values
    DERIVED_OPTION
    // counter keys
    N_RECORDS
    N_BY_KEY
    N_BY_VALUE
);
const MY_CONSTANT: u8 = 1; // additional fixed values not exposed as options

// main function called by xxx_tools main()
pub fn stream() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[MY_OPTION]);
    cfg.set_bool_env(&[MY_FLAG]);

    // set derived config values
    cfg.set_usize(DERIVED_OPTION, if cfg.equals_string(MY_OPTION, "some_value") { 22 } else { 33 });

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_RECORDS, "number of records processed"), // add regular counters, as many as needed
    ]);
    // ctrs.add_keyed_counters(&[
    //     (N_BY_KEY, "count of unique keys") // add keyed counters, as many as needed, where each counter counts categorical keys
    // ]);
    // ctrs.add_indexed_counters(&[
    //     (N_BY_VALUE, "distribution of integer values") // add indexed counters, as many as needed, where each counter counts integer indices
    // ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.print("initializing");

    // build tool support resources, i.e. data workers
    let mut tool = Tool {
        my_worker: 0,
        // continue with as many tool workers as needed
    };

    // process records in a stream
    w.log.print("doing something on streamed records");
    let mut rs = RecordStreamer::new();
    rs // choose whatever record_streamer and options match your input data format
        .group_by_in_place_serial(|alns: &mut [SamRecord]| record_parser(
            alns, 
            &mut w,
            &mut tool,
        ), &["grouping_field"]);

    // report counter values
    w.ctrs.print_all();
    Ok(())
}

// process a record or group of records
fn record_parser(
    records: &mut [SamRecord], 
    w: &mut Workflow,
    tool: &mut Tool,
) -> Result<Vec<usize>, Box<dyn Error>> {
    w.ctrs.increment(N_RECORDS);
    // w.ctrs.increment_keyed(N_BY_KEY, records[0].key_field.clone());

    // use w.cfg and tool workers to process each record  as needed

    // return the value appropriate for the streamer to indicate which records to print to STDOUT
    Ok((0..records.len()).collect())
}

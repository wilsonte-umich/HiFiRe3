//! Command-line tools for processing HiFiRe3 data.

// dependencies
use std::env;
use std::error::Error;

// modules
mod tools;
mod formats;
mod junctions;

// constants
const TOOLS_NAME: &str = "hf3_tools";

// load and process data
fn main() -> Result<(), Box<dyn Error>> {

    // read command line arguments
    let mut args: Vec<String> = env::args().collect();
    args = args[1..].to_vec(); // drop executable name 
    if args.len() == 0 { // check for something to do, i.e., a tool to run
        eprintln!("{}: missing tool or command", TOOLS_NAME);
        Err(format!("usage: {} <tool> [additional arguments]", TOOLS_NAME))?
    }
    let tool = args[0].clone(); // drop tool name
    // args = args[1..].to_vec();

    // dispatch to tool or command
    match tool.as_str() {

        /*--------------------------------------------------------------
        basecall ONT
        ------------------------------------------------------------- */
        // RE site-aware trimming for ONT reads
        "trim_ont" => tools::trim_ont::stream(),

        /*--------------------------------------------------------------
        analyze fragments
        ------------------------------------------------------------- */
        // sort and examine reads in the output stream from minimap2
        "analyze_alignments" => tools::analyze_alignments::stream(),

        // extract RE site positions from read 5' (and 3') ends
        "extract_endpoints" => tools::extract_endpoints::stream(),

        // // // create binary lookup files from extracted RE site positions
        // // "create_index" => tools::create_index::stream(),

        // // apply various alignment and SV quality filters and error correction mechanisms to create SITE_SAM
        // "apply_filters" => tools::apply_filters::stream(),

        /*--------------------------------------------------------------
        unrecognized pipeline action tool
        ------------------------------------------------------------- */
        _ => Err(format!("{}: unknown tool or command: {}", TOOLS_NAME, tool))?
    }
}

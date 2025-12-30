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

        // compression and reformatting of legacy ONT POD5 files
        "reformat_ont" => tools::reformat_ont::stream(),

        /*--------------------------------------------------------------
        basecall PacBio
        ------------------------------------------------------------- */
        // merge fwd and rev PacBio reads into pbFree "unleaded" reads
        "basecall_pacbio" => tools::basecall_pacbio::stream(),

        /*--------------------------------------------------------------
        analyze fragments
        ------------------------------------------------------------- */
        // sort (per read) and examine read alignments in the output stream from minimap2
        "analyze_alignments" => tools::analyze_alignments::stream(),

        // assess inserts based on RE site matching and insert sizing/stem lengths
        "analyze_inserts" => tools::analyze_inserts::stream(),

        /*--------------------------------------------------------------
        analyze SVs
        ------------------------------------------------------------- */
        // analyze reads for evidence of shared SV junctions
        "analyze_svs" => tools::analyze_svs::stream(),

        /*--------------------------------------------------------------
        unrecognized pipeline action tool
        ------------------------------------------------------------- */
        _ => Err(format!("{}: unknown tool or command: {}", TOOLS_NAME, tool))?
    }
}

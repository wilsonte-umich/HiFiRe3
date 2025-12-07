//! Crate `hf3_tools` provides various command-line tools for processing 
//! HiFiRe3 data. Most operate on standard input and output streams.

// dependencies
use std::env;
use std::error::Error;

// modules
mod tools;
mod formats;

// load and process data
fn main() -> Result<(), Box<dyn Error>> {

    // read command line arguments
    let args: Vec<String> = env::args().collect();
    if args.len() == 0 {
        Err("Usage: hf3_tools <command> [additional arguments]")?
    }

    // dispatch command
    match args[1].as_str() {

        // RE site-aware trimming for ONT reads
        "trim_ont" => tools::trim_ont::stream(),

        // "sv_apply_filters" => tools::sv_apply_filters::stream(),

        // unrecognized command
        _ => Err(format!("Unrecognized command: {}", args[1]))?
    }
}

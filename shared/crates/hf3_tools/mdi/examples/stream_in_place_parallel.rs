//! A simple app to show the use of mdi::stream::RecordStreamer::stream_in_place_parallel().
//! Compatible with output streamed from mdi_streamer/make_tsv.pl.

// dependencies
// use std::{thread, time};
// use rand::Rng;
use mdi::stream::RecordStreamer;
use serde::{Deserialize, Serialize};

// structures, with support for record parsing using serde
#[derive(Serialize, Deserialize)]
struct MyRecord {
    group:  u32,
    record: u32,
    name:   String,
    random: u32,
}

// constants, for parallel processing
const METHOD:      &str  = "stream_in_place_parallel";
const N_CPU:       usize = 4;
const BUFFER_SIZE: usize = 1000; // number of records to buffer before parallel processing

// main
fn main() {

    // demonstrate passing of immutable values to the record parser
    let proof: String = METHOD.to_string();
    let record_parser = |input_record: &mut MyRecord| -> Option<()> {
        parse_with_proof(input_record, &proof)
    };
    RecordStreamer::new()
        .stream_in_place_parallel(record_parser, N_CPU, BUFFER_SIZE);
}

// record parsing function
// records are updated by reference, returning None or Some(()) to enact filtering
fn parse_with_proof(input_record: &mut MyRecord, proof: &str) -> Option<()> {

    // // simulate a slow process by sleeping for a random number of milliseconds
    // // output order will be retained by par_iter.map()
    // let milli_seconds: u64 = rand::thread_rng().gen_range(0..5);
    // thread::sleep(time::Duration::from_millis(milli_seconds)); 

    // filter against some records by returning None
    if input_record.group > 5 && input_record.group < 10 {
        None

    // update the remaining records to show we did something
    } else {
        input_record.random *= 100;
        input_record.name = format!("{}-{}", input_record.name, proof);

        // return Some(()) to indicate success
        // do not need to return the record since it is updated in place
        Some(())
    }
}

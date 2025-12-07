//! A simple app to show the use of mdi::stream::RecordStreamer::group_by_in_place_parallel().
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
const METHOD:      &str  = "group_by_in_place_parallel";
const N_CPU:       usize = 4;
const BUFFER_SIZE: usize = 1000;

// main
fn main() {

    // demonstrate passing of immutable values to the record parser
    let proof: String = METHOD.to_string();
    let record_parser = |input_record: &mut Vec<MyRecord>| -> Option<()> {
        parse_with_proof(input_record, &proof)
    };

    // in this example, we group by a single field
    RecordStreamer::new()
        .group_by_in_place_parallel(record_parser, &["group"], N_CPU, BUFFER_SIZE);
}

// record parsing function
// records are updated by reference, returning None or Some(()) to enact filtering at the group level
fn parse_with_proof(input_record_group: &mut Vec<MyRecord>, proof: &str) -> Option<()> {

    // // simulate a slow process by sleeping for a random number of milliseconds
    // // output order will be retained by par_iter.map()
    // let milli_seconds: u64 = rand::thread_rng().gen_range(0..5);
    // thread::sleep(time::Duration::from_millis(milli_seconds)); 

    // filter against some record groups by returning None
    let group = input_record_group[0].group;
    if group > 5 && group < 10 {
        None

    // update the remaining records to show we did something
    } else {
        for input_record in input_record_group.iter_mut() {
            input_record.random *= 100;
            input_record.name = format!("{}-{}-{}", input_record.name, proof, group);
        }

        // return Some(()) to indicate success
        // do not need to return the record since it is updated in place
        Some(())
    }
}

//! A simple app to show the use of mdi::stream::RecordStreamer::stream_in_place_serial().
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

// main
fn main() {
    RecordStreamer::new()
        .stream_in_place_serial(record_parser);
}

// record parsing function
// records are updated by reference, returning None or Some(()) to enact filtering
fn record_parser(input_record: &mut MyRecord) -> Option<()> {

    // // simulate a slow process by sleeping for a random number of milliseconds
    // // output order will be retained (obligatorily since records are processed serially)
    // let milli_seconds: u64 = rand::thread_rng().gen_range(0..5);
    // thread::sleep(time::Duration::from_millis(milli_seconds)); 

    // filter against some records by returning None
    if input_record.group > 5 && input_record.group < 10 {
        None

    // update the remaining records to show we did something
    } else {
        input_record.random *= 100;
        input_record.name = format!("{}-{}", input_record.name, "stream_in_place_serial");

        // return Some(()) to indicate success
        // do not need to return the record since it is updated in place
        Some(())
    }
}

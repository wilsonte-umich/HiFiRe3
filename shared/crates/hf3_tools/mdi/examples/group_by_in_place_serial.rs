//! A simple app to show the use of mdi::stream::RecordStreamer::group_by_in_place_serial().
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
    // in this example, we group by a single field
    RecordStreamer::new()
        .group_by_in_place_serial(record_parser, &["group"]);
}

// record parsing function
// records are updated by reference, returning None or Some(()) to enact filtering at the group level
fn record_parser(input_record_group: &mut Vec<MyRecord>) -> Option<()> {

    // // simulate a slow process by sleeping for a random number of milliseconds
    // // output order will be retained (obligatorily since records are processed serially)
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
            input_record.name = format!("{}-{}-{}", input_record.name, "group_by_in_place_serial", group);
        }

        // return Some(()) to indicate success
        // do not need to return the record since it is updated in place
        Some(())
    }
}

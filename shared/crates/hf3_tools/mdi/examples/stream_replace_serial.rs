//! A simple app to show the use of mdi::stream::RecordStreamer::stream_replace_serial().
//! Compatible with output streamed from mdi_streamer/make_tsv.pl.

// dependencies
// use std::{thread, time};
// use rand::Rng;
use mdi::stream::RecordStreamer;
use serde::{Deserialize, Serialize};

// structures, with support for record parsing using serde
#[derive(Serialize, Deserialize)]
struct InputRecord {
    group:  u32,
    record: u32,
    name:   String,
    random: u32,
}
#[derive(Serialize, Deserialize)]
struct OutputRecord {
    group:  u32,
    record: u32,
    name:   String,
    random: u32,
    proof:  String, // to illustrate adding a field to the output record
}

// implement transformation (cast) from input to output records
// this is the slowest step as it requires copying the input data to a new output record
impl From<&InputRecord> for OutputRecord {
    fn from(input_record: &InputRecord) -> Self {
        OutputRecord {
            group:  input_record.group,        // integers implement Copy so no clone is required
            record: input_record.record,
            name:   input_record.name.clone(), // clone is required since name is a String
            random: input_record.random,
            proof:  "".to_string(),            // initialise proof to empty string
        }
    }
}

// main
fn main() {
    RecordStreamer::new()
        .stream_replace_serial(record_parser);
}

// record parsing function
// input records are immutable and must be transformed to output records
fn record_parser(input_record: &InputRecord) -> Option<Vec<OutputRecord>> {

    // // simulate a slow process by sleeping for a random number of milliseconds
    // // output order will be retained (obligatorily since records are processed serially)
    // let milli_seconds: u64 = rand::thread_rng().gen_range(0..5);
    // thread::sleep(time::Duration::from_millis(milli_seconds)); 

    // filter against some records by returning None
    if input_record.group > 5 && input_record.group < 10 {
        None
    } else {

        // create a new output record
        let mut output_record = OutputRecord::from(input_record);
        // let mut output_record: OutputRecord = input_record.into(); // alternative syntax

        // update the proof fields
        output_record.random *= 100;
        output_record.proof = format!("{}-{}", output_record.name, "stream_replace_serial");

        // return the new output record(s)
        // returning a vector of records transfers metadata ownership to RecordStreamer
        // without a deep copy of the allocated record data on the heap
        Some(vec![output_record])
    }
}

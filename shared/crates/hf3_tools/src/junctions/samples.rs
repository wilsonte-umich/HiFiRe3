//! Handle sample_bit to sample_name conversion, especially during merging.

// dependencies
use std::error::Error;
use serde::{Serialize, Deserialize};
use mdi::pub_key_constants;
use mdi::workflow::Config;
use mdi::{InputCsv, OutputCsv};
use super::junction::FinalJunction;

// constants
pub_key_constants!(
    // from environment variables
    SV_SAMPLES_FILE
);

/// Relate sample bit to name.
#[derive(Serialize, Deserialize)]
pub struct Sample {
    sample_bit:  u32,
    sample_name: String,
}
impl Sample {
    /// Update sample bit from multiple merging files.
    pub fn parse_merge_samples(
        cfg: &mut Config, 
        merge_input_dirs: &[&str]
    ) -> Result<Vec<FinalJunction>, Box<dyn Error>>{
        let mut samples_file = OutputCsv::open_env(cfg, SV_SAMPLES_FILE, None); 
        let mut sample_slots_used_all: u32 = 0;
        let mut jxns: Vec<FinalJunction> = Vec::new();
        for dir in merge_input_dirs {
            eprintln!("  {}", dir);

            // input and re-output sample from one merge source
            let mut sample_slots_used_dir = 0;
            let mut reader = InputCsv::open_file_from_glob(
                dir, "sv_samples.txt", 
                b'\t', true
            )?;
            for result in reader.deserialize(){
                let mut sample: Sample = result?;
                sample.sample_bit = sample.sample_bit << sample_slots_used_all;
                samples_file.serialize(&sample);
                sample_slots_used_dir += 1;
            }

            // bit shift the sample_bits in the recorded jxn to match above
            // bit shift applies the same to all source samples, whose order is maintained above
            let mut reader = InputCsv::open_file_from_glob(
                dir, "final_junctions_1.txt.bgz", 
                b'\t', true
            )?;
            for result in reader.deserialize(){
                let mut jxn: FinalJunction = result?;
                jxn.sample_bits = jxn.sample_bits << sample_slots_used_all;
                jxns.push(jxn);
            }

            // prepare for the next merge source
            sample_slots_used_all += sample_slots_used_dir;  
        }  

        // return
        Ok(jxns)    
    }
}

//! Handle the input BAM stream; dispatch strand pairs for merging.

// modules
mod usable;

// dependencies
use std::error::Error;
use std::path::Path;
use std::collections::HashMap;
use std::str::from_utf8_unchecked;
use crossbeam::channel::Sender;
use rust_htslib::{bam::{Read, Reader, Record as BamRecord}};
use mdi::pub_key_constants;
use mdi::workflow::Config;
use crate::formats::hf_tags::*;
use super::{BufferedStrand, StrandPair, MergeResult, bam};

/// StrandBuffer structure for caching PacBio strands
/// while waiting for their matching strand.
struct StrandBuffer {
    movies: Vec<String>,
    pub strand_buffer: HashMap<usize, BufferedStrand>,  
}
impl StrandBuffer {
    /// Create a new, empty StrandBuffer structure.
    pub fn new() -> Self {
        StrandBuffer {
            movies: Vec::new(),
            strand_buffer: HashMap::new(),
        }
    }
    /// Convert a movie name and ZMW hole number into a unique usize strand buffer key.
    pub fn get_key(&mut self, movie: &str, zmw: usize) -> usize {
        let movie_idx = if let Some(idx) = self.movies.iter().position(|m| m == movie) {
            idx
        } else {
            self.movies.push(movie.to_string());
            self.movies.len() - 1
        };
        movie_idx << 32 | zmw
    }
    /// Push a first encountered strand into the buffer for later merging.
    pub fn insert(&mut self, key: usize, usable: bool, ff:u8, ec: f32, strand: &BamRecord) {
        self.strand_buffer.insert(key, BufferedStrand {
            usable,
            ff,
            ec,
            seq:      strand.seq().as_bytes().iter().map(|&c| c as char).collect(),
            qual:     strand.qual().to_vec(),
            ip:       bam::get_tag_u8_vec_opt(strand, INTER_PULSE_DURATION),
            pw:       bam::get_tag_u8_vec_opt(strand, PULSE_WIDTH),
        });
    }
    /// Attempt to remove a first strand from the buffer for merging to a second strand.
    pub fn remove(&mut self, key: &usize) -> Option<BufferedStrand> {
        self.strand_buffer.remove(&key)
    }
}

// constants
pub_key_constants!(
    // from environment variables
    FAIL_BAM_FILE
    HIFI_BAM_FILE
    // read outcomes
    UNUSABLE_READ
    ONE_USABLE_STRAND
);

// called by crossbeam scope to create the single input thread
pub fn stream_bam_files(
    tx_strand_pair:  Sender<StrandPair>,
    tx_merge_result: Sender<MergeResult>,
) -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[FAIL_BAM_FILE, HIFI_BAM_FILE]);

    // instantiate buffering
    let mut this_strand = BamRecord::new();
    let mut strand_buffer = StrandBuffer::new();

    // scan the fail and hifi BAM files on order to build a cache of by-strand reads
    // merge strands when both are encountered and usable

    // for bam_file_key in &[FAIL_BAM_FILE, HIFI_BAM_FILE] {
    for bam_file_key in &[HIFI_BAM_FILE] {

        let bam_file = cfg.get_string(bam_file_key);
        let bam_file_name = Path::new(&bam_file).file_name().unwrap().to_str().unwrap();
        eprintln!("processing {}", bam_file_name);
        let mut bam = Reader::from_path(&bam_file)?;
        while let Some(r) = bam.read(&mut this_strand) {
            r.expect("Failed to parse BAM record");
            parse_record(
                &this_strand, &mut strand_buffer, 
                &tx_strand_pair, &tx_merge_result
            )?;
        }
    }
    Ok(())
}

// dispatch one strand of a PacBio by-strand read for processing
fn parse_record(
    this_strand:     &BamRecord,
    strand_buffer:   &mut StrandBuffer,
    tx_strand_pair:  &Sender<StrandPair>,
    tx_merge_result: &Sender<MergeResult>
) -> Result<(), Box<dyn Error>> {

    // parse ZMW from read name: m64011_190714_120746/14/ccs/rev
    // fail if the QNAME or hole number can't be parsed
    let qname = unsafe { from_utf8_unchecked(this_strand.qname()) };
    let parts: Vec<&str> = qname.split('/').collect();
    if parts.len() < 4 { 
        return Err(format!("Unexpected read name format: {}", qname).into());
    }
    let zmw = parts[1].parse::<u32>().map_err(|e|
        format!("Failed to parse ZMW from read name {}: {}", qname, e)
    )?;

    // assemble a buffer key for this strand's ZMW that is unique across movies
    let buffer_key = strand_buffer.get_key(parts[0], zmw as usize);

    // assess whether this strand is usable for merging
    let (this_usable, this_ff, this_ec) = usable::is_usable(this_strand);

    // process the second strand of a matching by-strand read pair
    // remove the cached data from the HashMap for memory management
    if let Some(prev_strand) = strand_buffer.remove(&buffer_key){
        let qname = format!("{}/{}/ccs", parts[0], zmw).into_bytes();
        let ff = prev_strand.ff | this_ff; // read reports both strands' fail flags

        // merge two usable strands into a single read
        // used downstream for both SNV and SV detection
        if prev_strand.usable && this_usable {
            tx_strand_pair.send(StrandPair{
                qname,
                ff,
                this: BufferedStrand {
                    usable: this_usable,
                    ff,
                    ec:     this_ec,
                    seq:    this_strand.seq().as_bytes().iter().map(|&c| c as char).collect(),
                    qual:   this_strand.qual().to_vec(),
                    ip:     bam::get_tag_u8_vec_opt(this_strand, INTER_PULSE_DURATION),
                    pw:     bam::get_tag_u8_vec_opt(this_strand, PULSE_WIDTH),
                },
                prev: prev_strand,
            })?;

        // transmit reads with one usable strand directly to writer without the two-strand tag
        // used downstream for SV detection only
        } else if prev_strand.usable {
            tx_merge_result.send(MergeResult{
                outcome: ONE_USABLE_STRAND,
                qname,
                seq:  prev_strand.seq,
                qual: prev_strand.qual,
                ff,
                ec:   prev_strand.ec,
                dd:   None,
                sk:   None,
                dt:   None,
            })?;
        } else if this_usable {
            let seq: String = this_strand.seq().as_bytes().iter().map(|&c| c as char).collect();
            tx_merge_result.send(MergeResult{
                outcome: ONE_USABLE_STRAND,
                qname,
                seq,
                qual: this_strand.qual().to_vec(),
                ff,
                ec: this_ec,
                dd: None,
                sk: None,
                dt: None,
            })?;

        // reject reads with no usable strands due to adapter parsing failures, low coverage, etc.
        } else {
            tx_merge_result.send(MergeResult{
                outcome: UNUSABLE_READ,
                qname,
                seq: String::new(),
                qual: Vec::new(),
                ff,
                ec: 0.0,
                dd: None,
                sk: None,
                dt: None,
            })?;
        }

    // cache the needed bits of the first encountered strand of a ZMW while waiting for the second
    } else {
        strand_buffer.insert(buffer_key, this_usable, this_ff, this_ec, this_strand);
    }
    Ok(())
}


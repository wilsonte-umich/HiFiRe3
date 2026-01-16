//! Merge the two strands of PacBio CCS reads into a single sequence.
// input:
//     "fail" and "hifi" CCS reads as FAIL_BAM_FILE, HIFI_BAM_FILE
//      called with the by-strand option, and optionally with kinetics tags
// output:
//     READ_BAM_FILE, with:
//         - one read per original ZMW, either:
//              - merging both strands for error corrected SNV and SV detection
//              - passing a single usable strand as is for SV detection only
//              - suppressing reads with no usable strands
//     PACBIO_BASECALL_KINETICS file with kinetics values stratified by 3-base context

// modules
mod channel_reader;
mod channel_worker;
mod channel_kinetics;

// dependencies
use std::error::Error;
use crossbeam::channel::bounded;
use minimap2::{Aligner as Minimap2};
use faimm::IndexedFasta;
use rust_htslib::bam::{Reader, Read, Writer, Record as BamRecord, record::Aux, Header, Format};
use mdi::pub_key_constants;
use mdi::workflow::{Config, Counters};
use crate::formats::hf3_tags::*;
use channel_kinetics::kinetics::{StrandMetadata, FrameCounts};

/// BufferedStrand holds the data needed for a single PacBio strand.
struct BufferedStrand {
    pub usable:   bool,
    pub ff:       u8, // same names as tags
    pub ec:       f32,
    pub seq:      String,
    pub qual:     Vec<u8>, 
    pub ip:       Option<Vec<u8>>, // inter-pulse durations
    pub pw:       Option<Vec<u8>>, // pulse widths
}

// MergePair holds the data needed to merge two strands.
struct StrandPair {
    qname: Vec<u8>,
    ff:    u8,
    this:  BufferedStrand,
    prev:  BufferedStrand,
}

/// KineticsInstance supports kinetics data collection in the side-effect channel.
struct KineticsInstance {
    strand_metadata: StrandMetadata,
    frame_counts:    FrameCounts,
}

/// MergeResult carries the outcome of attempting to merge two strands.
struct MergeResult {
    outcome: &'static str, 
    qname:   Vec<u8>,
    seq:     String,
    qual:    Vec<u8>,
    ff:      u8,
    ec:      f32,
    dd:      Option<String>,
    sk:      Option<Vec<u16>>,
    dt:      Option<u8>,
}

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "basecall_pacbio";
pub_key_constants!(
    // from environment variables
    N_CPU
    HIFI_BAM_FILE // here for header copying
    READ_BAM_FILE // the output BAM file for merged reads
    MINIMAP2_INDEX_WRK
    GENOME_FASTA
    // counter keys
    N_READS
    N_READS_BY_OUTCOME
    N_STRAND_DIFFERENCE_TYPES
    N_READS_BY_MINIMAP2
    // merge outcomes
    UNUSABLE_READ
);
const CHANNEL_CAPACITY: usize = 16; // channel buffer size
const MIN_READ_LEN: usize = 250;
const MAX_READ_LEN: usize = 25000;

// main basecall pacbio function called by hf3_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_usize_env(&[N_CPU]);
    cfg.set_string_env(&[HIFI_BAM_FILE, READ_BAM_FILE, GENOME_FASTA, MINIMAP2_INDEX_WRK]);

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_READS, "number of reads with two opposing strands"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_READS_BY_OUTCOME, "number of reads by strand merging outcome"),
        (N_STRAND_DIFFERENCE_TYPES, "number of merged reads by strand difference type"),
        (N_READS_BY_MINIMAP2, "number of reads by minimap2 strand alignment outcome"),
    ]);

    // instantiate shared minimap2 aligner
    let mm2_index = cfg.get_string(MINIMAP2_INDEX_WRK).to_string();
    let minimap2 = Minimap2::builder()
        .map_hifi()
        .with_index(mm2_index, None)
        // .with_cigar()
        .expect("Failed to build minimap2 index");

    // instantiate the FAI sequence fetcher
    let genome_fasta = cfg.get_string(GENOME_FASTA).to_string();
    let fa = IndexedFasta::from_file(&genome_fasta)
        .expect("Error opening genome FASTA file");

    // read the header from HIFI_BAM_FILE to use for the output
    let hifi_path = cfg.get_string(HIFI_BAM_FILE).to_string();
    let hifi_bam = Reader::from_path(&hifi_path)?;
    let header = Header::from_template(hifi_bam.header());

    // create output BAM header and writer
    let output_path = cfg.get_string(READ_BAM_FILE).to_string();
    let mut bam_writer = Writer::from_path(&output_path, &header, Format::Bam)?;

    // create channels for fan-in/fan-out parallel processing, and kinetics side-effect
    let (tx_strand_pair, rx_strand_pair)    
        = bounded::<StrandPair>(CHANNEL_CAPACITY);
    let (tx_kinetics, rx_kinetics) 
        = bounded::<KineticsInstance>(CHANNEL_CAPACITY);
    let (tx_merge_result, rx_merge_result) 
        = bounded::<MergeResult>(CHANNEL_CAPACITY);

    // spawn threads
    crossbeam::scope(|scope| {

        // reader: read strands from input BAMs and dispatch for merging
        let tx_mr  = tx_merge_result.clone();
        scope.spawn(move |_| {
            channel_reader::stream_bam_files(
                tx_strand_pair, 
                tx_mr, // for short-circuting merge worker
            ).unwrap();
        });

        // workers: merge strand pairs from reader in parallel
        // create various types of merge results and kinetics side-effects
        for _ in 0..*cfg.get_usize(N_CPU) {
            let rx_strand_pair  = rx_strand_pair.clone();
            let tx_kinetics = tx_kinetics.clone();
            let tx_merge_result  = tx_merge_result.clone();
            scope.spawn(|_| {
                channel_worker::merge_strand_pairs(
                    rx_strand_pair, 
                    tx_kinetics,
                    tx_merge_result,
                    &minimap2,
                    &fa,
                ).unwrap();
            });
        }

        // kinetics side-effect receiver channel
        scope.spawn(move |_| {
            channel_kinetics::collect_kinetics(
                rx_kinetics
            ).unwrap();
        });

        // clean up channels
        drop(tx_kinetics);
        drop(tx_merge_result); 

        // process merge results as they arrive
        for merge_result in rx_merge_result {

            // update counters
            ctrs.increment(N_READS);
            ctrs.increment_keyed(N_READS_BY_OUTCOME, merge_result.outcome);

            // write merged read if at least one usable strand
            if merge_result.outcome != UNUSABLE_READ {
                let mut read = BamRecord::new();
                read.set(
                    &merge_result.qname, 
                    None, 
                    merge_result.seq.as_bytes(), 
                    &merge_result.qual
                );
                read.push_aux(
                    PACBIO_FAIL.as_bytes(), 
                    Aux::U8(merge_result.ff)
                ).unwrap();
                read.push_aux(
                    PACBIO_EFF_COVERAGE.as_bytes(), 
                    Aux::Float(merge_result.ec)
                ).unwrap();
                if let Some(dd) = merge_result.dd {  
                    read.push_aux(
                        STRAND_DIFFERENCES.as_bytes(), 
                        Aux::String(&dd)
                    ).unwrap();
                }
                if let Some(sk) = merge_result.sk { 
                    read.push_aux(
                        SUBSTITUTION_KINETICS.as_bytes(), 
                        Aux::ArrayU16(sk.as_slice().into())
                    ).unwrap();
                }
                if let Some(dt) = merge_result.dt {
                    read.push_aux(
                        STRAND_DIFFERENCE_TYPES.as_bytes(),
                        Aux::U8(dt)
                    ).unwrap();
                    ctrs.increment_keyed(
                        N_STRAND_DIFFERENCE_TYPES, 
                        dt.to_string().as_str()
                    );
                }
                bam_writer.write(&read).unwrap();
            }
        }
    }).expect("Crossbeam scope panicked");

    // print counts
    ctrs.print_grouped(&[
        &[N_READS],
        &[N_READS_BY_OUTCOME],
        &[N_STRAND_DIFFERENCE_TYPES],
    ]);
    Ok(())
}

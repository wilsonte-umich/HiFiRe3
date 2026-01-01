//! Merge the two strands of PacBio CCS reads into a single sequence.
// input:
//     "fail" and "hifi" CCS reads as FAIL_BAM_FILE, HIFI_BAM_FILE
//      called with the by-strand option, and optionallys with kinetics tags
// output:
//     READ_BAM_FILE, with:
//         - one read per original ZMW, either:
//              - merging both strands for error corrected SNV and SV detection
//              - passing a single usable strand as is for SV detection only
//              - suppressing reads with no usable strands
//     PACBIO_BASECALL_KINETICS file with kinetics values stratified by 3-base context

// modules
mod bam;
mod usable;
mod strand_merger;
mod kinetics;

// dependencies
use std::error::Error;
use std::str::from_utf8_unchecked;
use rust_htslib::bam::{Read, Reader, Writer, Record as BamRecord, record::Aux, Header, Format};
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use genomex::sequence::{rc_acgt_str, Aligner, AlignmentStatus};
use crate::formats::hf_tags::*;
use strand_merger::{BufferedStrand, StrandBuffer, StrandMerger};

// Tool structure, for passing data workers to record processing functions.
struct Tool {
    strand_buffer: StrandBuffer,
    aligner:       Aligner,
    strand_merger: StrandMerger,
    bam_writer:    Writer,
}

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "basecall_pacbio";
pub_key_constants!(
    // from environment variables
    FAIL_BAM_FILE
    HIFI_BAM_FILE
    READ_BAM_FILE // the output BAM file for merged reads
    PACBIO_BASECALL_KINETICS // output kinetics file path
    // counter keys
    N_STRANDS
    N_READS
    N_READS_BY_OUTCOME
    // read outcomes
    UNUSABLE_READ
    ONE_USABLE_STRAND
    TWO_STRANDS_MERGED
    FAILED_STRAND_ALIGNMENT
);
const MAX_SHIFT: usize = 10; // allow only small total indel shifts when aligning strands

// main basecall pacbio function called by hf3_tools main()
pub fn stream() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[FAIL_BAM_FILE, HIFI_BAM_FILE, READ_BAM_FILE, PACBIO_BASECALL_KINETICS]);

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_STRANDS, "number of strand records processed"), 
        (N_READS,   "number of reads with two opposing strands"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_READS_BY_OUTCOME, "number of reads by outcome"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.print("initializing");

    // read the header from HIFI_BAM_FILE to use for the output
    let hifi_path = w.cfg.get_string(HIFI_BAM_FILE).to_string();
    let hifi_bam = Reader::from_path(&hifi_path)?;
    let header = Header::from_template(hifi_bam.header());

    // create output BAM header and writer
    let output_path = w.cfg.get_string(READ_BAM_FILE).to_string();
    let bam_writer = Writer::from_path(&output_path, &header, Format::Bam)?;

    // build tool support resources
    let mut tool = Tool {
        strand_buffer: StrandBuffer::new(),
        aligner: Aligner::new(
            usable::MAX_READ_LEN + 1, 
            usable::MAX_READ_LEN + 1
        ).max_shift(MAX_SHIFT),
        strand_merger: StrandMerger::new(),
        bam_writer,
    };

    // scan the fail and hifi BAM files to build a cache of by-strand reads
    // merge strands when both are encountered and usable
    let mut this_strand = BamRecord::new();
    stream_records(FAIL_BAM_FILE, &mut this_strand, &mut w, &mut tool)?;
    stream_records(HIFI_BAM_FILE, &mut this_strand, &mut w, &mut tool)?;

    // finalize kinetics file
    tool.strand_merger.kinetics.write_kinetics_file(
        w.cfg.get_string(PACBIO_BASECALL_KINETICS)
    )?;

    // print counts
    w.ctrs.print_grouped(&[
        &[N_STRANDS, N_READS],
        &[N_READS_BY_OUTCOME],
    ]);
    Ok(())
}

// stream records from fail then hifi unmapped BAM files
fn stream_records(
    file_key:    &str,
    this_strand: &mut BamRecord,
    w:           &mut Workflow,
    tool:        &mut Tool,
) -> Result<(), Box<dyn Error>> {
    let file_path = w.cfg.get_string(file_key).to_string();
    let mut bam = Reader::from_path(&file_path)?;
    while let Some(r) = bam.read(this_strand) {
        r.expect("Failed to parse BAM record");
        parse_record(this_strand, w, tool)?;
    }
    Ok(())
}

// dispatch one strand of a PacBio by-strand read for processing
fn parse_record(
    this_strand: &mut BamRecord,
    w:           &mut Workflow,
    tool:        &mut Tool,
) -> Result<(), Box<dyn Error>> {
    w.ctrs.increment(N_STRANDS);

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
    let buffer_key = tool.strand_buffer.get_key(parts[0], zmw as usize);

    // assess whether this strand is usable for merging
    let (this_usable, this_ec) = usable::is_usable(this_strand);

    // process the second strand of a matching by-strand read pair
    // remove the cached data from the HashMap for memory management
    if let Some(prev_strand) = tool.strand_buffer.remove(&buffer_key){
        w.ctrs.increment(N_READS);
        let qname = format!("{}/{}/ccs", parts[0], zmw);
        let qname = qname.as_bytes(); // qname without /fwd or /rev suffix

        // merge two usable strands into a single read
        // used downstream for both SNV and SV detection
        if prev_strand.usable && this_usable {
            merge_strands(this_strand, prev_strand, qname, w, tool)?;
        
        // pass reads with one usable strand as is without the two-strand tag set
        // used downstream for SV detection only
        } else if prev_strand.usable {
            commit_read(
                ONE_USABLE_STRAND, qname, 
                &prev_strand.seq, &prev_strand.qual,
                prev_strand.ec, None, 
                w, &mut tool.bam_writer
            )?;
        } else if this_usable {
            let seq: String = this_strand.seq().as_bytes().iter().map(|&c| c as char).collect();
            commit_read(
                ONE_USABLE_STRAND, qname, 
                &seq, this_strand.qual(), 
                this_ec, None,
                w, &mut tool.bam_writer
            )?;

        // reject reads with no usable strands due to:
        //  - adapter parsing failures
        //  - low effective coverage
        } else {
            w.ctrs.increment_keyed(N_READS_BY_OUTCOME, UNUSABLE_READ);
        }

    // cache the essential parts of the first encountered strand 
    // for the ZMW while waiting for the second strand
    } else {
        tool.strand_buffer.insert(buffer_key, this_usable, this_ec, this_strand);
    }
    Ok(())
}

// merge two usable strands into a single read, masking disagreements to Ns
fn merge_strands(
    this_strand: &mut BamRecord,
    prev_strand: BufferedStrand,
    qname:       &[u8],
    w:           &mut Workflow,
    tool:        &mut Tool,
) -> Result<(), Box<dyn Error>> {
    let this_ec = bam::get_tag_f32(this_strand, PACBIO_EFF_COVERAGE);

    // reverse complement the previous strand sequence for alignment
    let qry = rc_acgt_str(&prev_strand.seq);
    let tgt = this_strand.seq().as_bytes().iter().map(|&c| c as char).collect::<String>();

    // align the two strands
    let aln = tool.aligner.align(&qry, &tgt, None, true);

    // handle the unxpected outcome of failed alignment between strands
    if aln.status != AlignmentStatus::AlignmentFound {
        return commit_read(
            FAILED_STRAND_ALIGNMENT, qname, 
            &tgt, this_strand.qual(),
            this_ec, None,
            w, &mut tool.bam_writer
        );
    }

    // merge the two strands into a single read
    let tgt_ip = bam::get_tag_u8_vec_opt(this_strand, INTER_PULSE_DURATION);
    let tgt_pw = bam::get_tag_u8_vec_opt(this_strand, PULSE_WIDTH);
    tool.strand_merger.merge_strands(
        &aln, &tgt, &prev_strand.seq, 
        tgt_ip, tgt_pw, prev_strand.ip, prev_strand.pw
    );

    // finalize the merged read with the two-strand tag set
    commit_read(
        TWO_STRANDS_MERGED, qname, 
        &tool.strand_merger.seq, &tool.strand_merger.qual,
        prev_strand.ec + this_ec, Some(&tool.strand_merger.dd_tag),
        w, &mut tool.bam_writer
    )
}

// finish processing and print a read with one reported single or merged strand
fn commit_read(
    outcome:    &str,
    qname:      &[u8],
    seq:        &str,
    qual:       &[u8],
    ec:         f32,
    dd:         Option<&str>,
    w:          &mut Workflow,
    bam_writer: &mut Writer,
) -> Result<(), Box<dyn Error>> {
    w.ctrs.increment_keyed(N_READS_BY_OUTCOME, outcome);
    let mut read = BamRecord::new();
    read.set(qname, None, seq.as_bytes(), qual);
    read.push_aux(PACBIO_EFF_COVERAGE.as_bytes(), Aux::Float(ec))?;
    if let Some(dd) = dd {  
        read.push_aux(STRAND_DIFFERENCES.as_bytes(), Aux::String(dd))?;
    }
    bam_writer.write(&read)?;
    Ok(())
}

//! Merge the two strands of PacBio CCS reads into a single sequence.
// input:
//     unaligned SAM stream on STDIN, specifically, PacBio CCS reads:
//         - called with the by-strand option
//         - including both "fail" and "hifi" CCS reads, provided in that order
// output:
//     updated SAM stream on STDOUT, with:
//         - one read per original ZMW, either:
//              - merging both strands for error corrected SNV and SV detection
//              - passing a single usable strand as is for SV detection only
//              - suppressing reads with no usable strands

// modules
mod bam;
mod usable;
mod strand_buffer;

// dependencies
use std::error::Error;
use std::collections::HashMap;
use std::str::from_utf8_unchecked;
use rust_htslib::bam::{Read, Reader, Record as BamRecord, record::Aux};
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use genomex::sequence::{rc_acgt_str, Aligner, AlignmentStatus};
use crate::formats::hf_tags::*;
use strand_buffer::BufferedStrand;

// Tool structure, for passing data workers to record processing functions.
struct Tool {
    strand_buffer: HashMap<u32, BufferedStrand>,
    aligner:       Aligner,
}

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "basecall_pacbio";
pub_key_constants!(
    // from environment variables
    FAIL_BAM_FILE
    HIFI_BAM_FILE
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
    cfg.set_string_env(&[FAIL_BAM_FILE, HIFI_BAM_FILE]);

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

    // build tool support resources
    let mut tool = Tool {
        strand_buffer: HashMap::new(),
        aligner: Aligner::new(
            usable::MAX_READ_LEN + 1, 
            usable::MAX_READ_LEN + 1
        ).max_shift(MAX_SHIFT),
    };

    // scan the fail and hifi BAM files to build a cache of by-strand reads
    // merge strands when both are encountered and usable
    let mut this_strand = BamRecord::new();
    stream_records(FAIL_BAM_FILE, &mut this_strand, &mut w, &mut tool)?;
    stream_records(HIFI_BAM_FILE, &mut this_strand, &mut w, &mut tool)?;

    // print counts
    // TODO
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

    // assess whether this strand is usable for merging
    let (this_usable, this_ec) = usable::is_usable(this_strand);

    // process the second strand of a matching by-strand read pair
    // remove the cached data from the HashMap for memory management
    if let Some(prev_strand) = tool.strand_buffer.remove(&zmw){
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
                prev_strand.seq, &prev_strand.qual,
                prev_strand.ec, None, 
                w
            )?;
        } else if this_usable {
            let seq = this_strand.seq().as_bytes().iter().map(|&c| c as char).collect();
            commit_read(
                ONE_USABLE_STRAND, qname, 
                seq, this_strand.qual(), 
                this_ec, None,
                w
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
        tool.strand_buffer.insert(zmw, BufferedStrand {
            usable:   this_usable,
            ec:       this_ec,
            seq:      this_strand.seq().as_bytes().iter().map(|&c| c as char).collect(),
            qual:     this_strand.qual().to_vec(),
            ip:       bam::get_tag_u8_vec(this_strand, INTER_PULSE_DURATION),
            pw:       bam::get_tag_u8_vec(this_strand, PULSE_WIDTH),
        });
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
            tgt, this_strand.qual(),
            this_ec, None,
            w
        );
    }

    // mask reference strand bases where the two strands disagree to Ns
    // create a tag roughly similar to the minimap2 cs tag to record the strand differences 
    // between-strand indels are always recorded as insertions to denote the possible 
    //     presence of the base (otherwise it would disappear from the output seq)
    let mut seq = String::new();        // the new merged sequence, with Ns at disagreement positions
    let mut dd_tag: String = String::new();     // the strand differences tag, recording bases at N positions in seq
    let mut ident_len: usize = 0;               // length of a current identical run; committed to dd_tag when a difference is encountered
    let mut del_buffer: String = String::new(); // buffer for deletion base series; committed to dd_tag when next match/insertion is encountered
    if aln.tgt_start0 > 0 { // maintain the same read length as the target strand for best RE-site matching
        seq.push_str('N'.to_string().repeat(aln.tgt_start0).as_str());
        dd_tag.push_str(&format!("?{}", aln.tgt_start0)); // ? operation for skipped bases at the start and end
    }
    aln.qry_on_tgt.iter().enumerate().for_each(|(j, qry)| {
        let i = aln.tgt_start0 + j;
        let tgt = &tgt[i..=i];
        if qry == tgt { // base identity between strands
            if !del_buffer.is_empty() {
                dd_tag.push_str(&format!("+{}", del_buffer)); // yes, +, see above
                del_buffer.clear();
            }
            seq.push_str(tgt);
            ident_len += 1;
        } else {
            if ident_len > 0 {
                dd_tag.push_str(&format!(":{}", ident_len.to_string()));
                ident_len = 0;
            }
            if qry.len() > 1 { // insertion on qry relative to tgt
                if !del_buffer.is_empty() {
                    dd_tag.push_str(&format!("+{}", del_buffer)); // yes, +, see above
                    del_buffer.clear();
                }
                // keep the base that the insertion bases were recorded on as is, it is always a match
                seq.push_str('N'.to_string().repeat(qry.len() - 1).as_str());
                dd_tag.push_str(&format!("+{}", qry[0..qry.len() - 1].to_string()));
                seq.push_str(tgt);
                ident_len = 1;
            } else if qry == "-" { // deletion on qry relative to tgt
                seq.push('N'); // again, bases were present on one strand but not the other
                del_buffer.push_str(tgt);
            } else { // base substitution between strands
                if !del_buffer.is_empty() {
                    dd_tag.push_str(&format!("+{}", del_buffer)); // yes, +, see above
                    del_buffer.clear();
                }
                seq.push('N'); // substitution one strand relative to the other
                // TODO: * tag includes trinucleotide context, and kinetics
                dd_tag.push_str(&format!("*{}{}", tgt, qry));
            }
        }
    });
    if ident_len > 0 { dd_tag.push_str(&format!(":{}", ident_len)); }
    if !del_buffer.is_empty() { dd_tag.push_str(&format!("+{}", del_buffer)); }
    if aln.tgt_end0 < tgt.len() - 1 {
        let n = tgt.len() - 1 - aln.tgt_end0;
        seq.push_str('N'.to_string().repeat(n).as_str());
        dd_tag.push_str(format!("?{}", n).as_str());
    }

    // create a two-level QUAL string with 
    //     Phred 40 = "I" for bases that agree between strands
    //     Phred  0 = "!" for disagreeing bases
    let qual: Vec<u8> = seq.chars().map(|b| 
        if b == 'N' { 0 } else { 40 }
    ).collect();

    // finalize the merged read with the two-strand tag set
    commit_read(
        TWO_STRANDS_MERGED, qname, 
        seq, &qual,
        prev_strand.ec + this_ec, Some(dd_tag),
        w
    )
}

// finish processing and print a read with one reported single or merged strand
fn commit_read(
    outcome:    &str,
    qname:      &[u8],
    seq:        String,
    qual:       &[u8],
    ec:         f32,
    dd:         Option<String>,
    w:          &mut Workflow,
) -> Result<(), Box<dyn Error>> {
    w.ctrs.increment_keyed(N_READS_BY_OUTCOME, outcome);
    let mut read = BamRecord::new();
    read.set(qname, None, seq.as_bytes(), qual);
    read.push_aux(PACBIO_EFF_COVERAGE.as_bytes(), Aux::Float(ec))?;
    if let Some(dd) = dd {  
        read.push_aux(STRAND_DIFFERENCES.as_bytes(), Aux::String(&dd))?;
    }
    // TODO: print to output BAM stream
    Ok(())
}

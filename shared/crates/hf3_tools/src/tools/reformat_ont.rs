//! action:
//!     for legacy file conversion:
//!         move trim metadata from QNAME to TAGS
//!     to reduce BAM file size:
//!         quantize ONT quality scores
//!         remove unneeded tags
//! input:
//!     unaligned SAM stream, specifically, a previously generated ONT UBAM file, on STDIN
//! output:
//!     updated SAM stream on STDOUT, with:
//!         - SEQ and QUAL length unchanged, i.e., no trimming performed by this tool
//!         - all prior tags removed, except retaining ML, MM, ch, rn, fn
//!         - new tag added: tl:Z: = adapter trim lengths in format `<5' trim>,<3' trim>`
//!         - QUAL field quantized to 8 levels to reduce BAM file size

// dependencies
use std::error::Error;
use mdi::RecordStreamer;
use genomex::sam::{SamRecord, SamQual};
use crate::formats::hf_tags::{TRIM_LENGTHS, StageTags};

// main ONT trim function called by hf3_tools main()
pub fn stream() -> Result<(), Box<dyn Error>> {

    // run ubam reformatting on each SAM record in a stream
    let mut rs = RecordStreamer::new();
    rs
        .comment(b'@')
        .no_trim()
        .flexible()
        .stream_in_place_serial(record_parser);

    // this tool runs silently and reports nothing
    Ok(())
}

// trim reads one at a time, on both ends
fn record_parser(aln: &mut SamRecord) -> Result<bool, Box<dyn Error>> {

    // quantize ONT quality scores to reduce BAM file size
    unsafe { SamQual::quantize_qual_scores(&mut aln.qual.qual); }

    // recover any current trim length metadata
    // first attempt to read from tl:Z: tag from QNAME
    // in legacy files, QNAME exited as QNAME:channel:trim5:trim3
    // e.g. 065a0f03-6ac4-4983-98f9-be92796eb0f9:800:35:0
    let qname_fields: Vec<&str> = aln.qname.as_str().split(':').collect();
    let (tl_tag_in, original_qname) = if qname_fields.len() == 4 {
        let trim5 = qname_fields[2];
        let trim3 = qname_fields[3];
        (
            Some(format!("{},{}", trim5, trim3)), 
            Some(qname_fields[0].to_string())
        )
    // otherwise try to read from tl:Z: tag
    // if None, we won't set the tl:Z: tag
    } else {
        (
            aln.get_tag_value(TRIM_LENGTHS), 
            None
        )
    };

    // remove unneeded tags to reduce BAM file size
    aln.tags.retain(&StageTags::BaseCalling.tag_added_by_stage()); 

    // add the tl:Z: tag indicating trim lengths, but only if found in input
    if let Some(tl_tag_in) = tl_tag_in {
        aln.tags.tags.push(format!("{}{}", TRIM_LENGTHS, tl_tag_in));
    }

    // restore the original QNAME
    if let Some(original_qname) = original_qname {
        aln.qname = original_qname;
    }

    // all input reads repeated to STDOUT
    Ok(true)
}

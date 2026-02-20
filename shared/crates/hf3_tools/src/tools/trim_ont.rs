//! action:
//!     find ONT adapters on previously untrimmed and unaligned read ends
//!     search includes blunt RE half-cut-site to enhance search power vs. adapter alone
//!     trim just the adapter portion, leaving the half-cut-site genomic bases intact
//!     quantize ONT quality scores to reduce BAM file size
//!     remove unneeded tags to reduce BAM file size
//!     add tag indicating trim status
//! input:
//!     unaligned SAM stream, specifically, a Dorado-compatible output stream, on STDIN
//! output:
//!     updated SAM stream on STDOUT, with:
//!         - SEQ and QUAL fields trimmed to remove adapters, when found
//!         - all prior tags removed, except retaining ML, MM, ch, rn, fn
//!         - new tag added: tl:Z: = adapter trim lengths in format `<5' trim>,<3' trim>`
//!         - QUAL field quantized to 8 levels to reduce BAM file size
//! notes:
//!     this script processes unaligned ONT reads, i.e., never reverse complemented, so first bases correspond to the 5' adapter
//!     all reads are expected to have 5' adapters, but some may be low quality bases, truncated, or chimeras split by Dorado
//!     3' adapters are always shorter, often truncated, and not always present for incompletely sequenced fragments
//!     there is only ever one SAM line per read, so no grouping is necessary

// dependencies
use std::error::Error;
use std::fs::read_to_string;
use mdi::pub_key_constants;
use mdi::workflow::Config;
use mdi::RecordStreamer;
use genomex::sam::{SamRecord, SamQual};
use genomex::sequence::{rc_acgt_str, Aligner, ForceQryTerminus};
use crate::formats::hf3_tags::{TRIM_LENGTHS, StageTags};

// tool support structures
struct Adapter {
    adapter:    String,
    adapter_re: String,
    exact:      String,
}

// constants
// const TOOL: &str = "trim_ont";
pub_key_constants!(
    // from environment variables
    BLUNT_RE_TABLE // adapter ligation HiFiRe3 only supports blunt REs
    ENZYME_NAME
    ADAPTER_SEQUENCE
    PLATFORM_MIN_INSERT_SIZE // reads shorter than this are not trimmed
    // derived values
    CUT_SITE_OFFSET
    CUT_SITE_LEFT_HALF
    CUT_SITE_RIGHT_HALF
    INCLUDE_RE_HALF_SITE
);
// minimum score for adapter detection
// 3 RE bases + 1 A-tail + 6 adapter bases
// also yields the number of bases to store for exact matching
const MIN_SCORE:       i32   = 10;
const SEARCH_SPACE_5:  usize = 60; // search spaces empirically determined from actual reads
const SEARCH_SPACE_3:  usize = 20; // wide enough to be sensitive, narrow enough to be fast and specific

// main ONT trim function called by hf3_tools main()
pub fn stream() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[BLUNT_RE_TABLE, ENZYME_NAME, ADAPTER_SEQUENCE]);
    cfg.set_usize_env(&[PLATFORM_MIN_INSERT_SIZE]);

    // build tool support resources
    initialize_restriction_enzyme(&mut cfg);
    let (adapter5, adapter3) = initialize_adapter_sequences(&mut cfg);
    let max_len = SEARCH_SPACE_5 * 2;
    let mut aligner = Aligner::new(
        max_len, 
        max_len
    ).suppress_alignment_map();

    // run ONT trimming on each SAM record in a stream
    RecordStreamer::new()
        .comment(b'@')
        .no_trim()
        .flexible()
        .stream_in_place_serial(|aln: &mut SamRecord| record_parser(
            aln, 
            &cfg, 
            &adapter5, 
            &adapter3, 
            &mut aligner,
        ));

    // because this tool is run repeatedly, it runs silently and reports nothing
    Ok(())
}

// initialize the RE site and associated config values
fn initialize_restriction_enzyme(cfg: &mut Config) {

    // initialize config values when there RE sites are not expected at read ends
    cfg.set_usize(CUT_SITE_OFFSET, 0);
    cfg.set_string_list(&[
        (CUT_SITE_LEFT_HALF,  "".to_string()),
        (CUT_SITE_RIGHT_HALF, "".to_string())
    ]);
    cfg.set_bool(INCLUDE_RE_HALF_SITE, cfg.get_string(ENZYME_NAME) != "NA");
    if !cfg.get_bool(INCLUDE_RE_HALF_SITE) { return; }

    // read blunt RE site metadata from table
    let blunt_re_table = cfg.get_string(BLUNT_RE_TABLE).to_string();
    let enzyme_name    = cfg.get_string(ENZYME_NAME).to_string();
    let content = read_to_string(&blunt_re_table)
        .unwrap_or_else(|e| panic!("could not open {}: {}", blunt_re_table, e));

    // loop through table to find specified enzyme
    for (i, line) in content.lines().enumerate() {
        if i == 0 { continue; } // skip header
        let fields: Vec<&str> = line.split(',').collect(); // enzyme,strand,cut_site,regex,offset,CpG_priority,...
        if fields[0] != enzyme_name { continue; }
        let cut_site = fields[2].to_uppercase();
        let cut_site_offset: usize = fields[4].parse::<usize>().unwrap_or_else(|e| 
            panic!("invalid offset in line {} in {}: {}", i + 1, blunt_re_table, e)
        );
        cfg.set_usize(CUT_SITE_OFFSET, cut_site_offset);
        cfg.set_string_list(&[
            (CUT_SITE_LEFT_HALF,  cut_site.chars().take(cut_site_offset).collect()),
            (CUT_SITE_RIGHT_HALF, cut_site.chars().skip(cut_site_offset).take(cut_site_offset).collect())
        ]);
    }

    // ensure that the enzyme was found if specified
    if cfg.get_string(CUT_SITE_LEFT_HALF).is_empty() {
        panic!("unrecognized enzyme {enzyme_name}; must be a blunt cutter in shared/modules/REs/blunt_enzymes.csv");
    }
}

// initalize the adapter sequences and associated config values
fn initialize_adapter_sequences(cfg: &mut Config) -> (Adapter, Adapter) {

    // duplex portion of the ONT kit adapter
    // for ligation kit, last T matches the one-base A-tail
    // fused to 5' genomic ends
    let adapter_core = cfg.get_string(ADAPTER_SEQUENCE).to_string(); 
    let mut adapter5 = Adapter {
        adapter:    adapter_core.clone(),
        adapter_re: format!("{}{}", adapter_core, cfg.get_string(CUT_SITE_RIGHT_HALF)),
        exact:      "".to_string(),
    };
    adapter5.exact = adapter5.adapter_re.chars()
        .skip(adapter5.adapter_re.len() - MIN_SCORE as usize)
        .collect();

    // reverse-complement adapter fused to 3' genomic ends
    let adapter_core_rc = rc_acgt_str(&adapter_core);
    let mut adapter3 = Adapter {
        adapter:    adapter_core_rc.clone(),
        adapter_re: format!("{}{}", cfg.get_string(CUT_SITE_LEFT_HALF), adapter_core_rc),
        exact:      "".to_string(),
    };
    adapter3.exact = adapter3.adapter_re.chars()
        .take(MIN_SCORE as usize)
        .collect();

    // return both adapters
    (adapter5, adapter3)
}

// trim reads one at a time, on both ends
fn record_parser(
    aln: &mut SamRecord, 
    cfg: &Config,
    adapter5: &Adapter,
    adapter3: &Adapter,
    aligner: &mut Aligner
) -> Result<bool, Box<dyn Error>> {

    // read too short to trim, filter out of final output
    let seq_len = aln.seq.len();
    if seq_len < *cfg.get_usize(PLATFORM_MIN_INSERT_SIZE) { return Ok(false); } 

    // execute adapter searches
    let (found5, trim_len5) = find_adapter(
        cfg,
        adapter5, 
        aligner,
        &aln.seq[..SEARCH_SPACE_5], 
        ForceQryTerminus::QryEnd
    );
    let (found3, trim_len3) = find_adapter(
        cfg,
        adapter3, 
        aligner,
        &aln.seq[seq_len - SEARCH_SPACE_3..], 
        ForceQryTerminus::QryStart
    );

    // trim reads to remove adapters
    if found3 { // order is important
        aln.seq       = aln.seq[      ..seq_len - trim_len3].to_string();
        aln.qual.qual = aln.qual.qual[..seq_len - trim_len3].to_string();
    }
    if found5 {
        aln.seq       = aln.seq[      trim_len5..].to_string();
        aln.qual.qual = aln.qual.qual[trim_len5..].to_string();
    }

    // quantize ONT quality scores to reduce BAM file size
    unsafe { SamQual::quantize_qual_scores(&mut aln.qual.qual); }

    // remove unneeded tags to reduce BAM file size
    aln.tags.retain(&StageTags::BaseCalling.tag_added_by_stage()); 

    // add the tl:Z: tag indicating trim lengths
    aln.tags.tags.push(format!(
        "{}{},{}", 
        TRIM_LENGTHS,
        if found5 { trim_len5 } else { 0 }, 
        if found3 { trim_len3 } else { 0 }
    ));

    // all input reads of sufficient length pass to output regardless of trimming status
    Ok(true)
}

// execute an efficient adapter search
// searches force alignment to the end of the adapter closest to the adapter-gDNA transition for specificity (also the highest quality bases)
fn find_adapter(
    cfg: &Config,
    adapter: &Adapter, 
    aligner: &mut Aligner,
    read_segment: &str, 
    forced_end: ForceQryTerminus, 
) -> (bool, usize) {
    let mut score: i32 = 0;
    let mut trim_len: usize = 0;
    let cut_site_offset = *cfg.get_usize(CUT_SITE_OFFSET);
    let include_re_half_site = *cfg.get_bool(INCLUDE_RE_HALF_SITE);

    // for speed, search first for the adapter-RE transition by exact matching
    if let Some(i0) = read_segment.find(&adapter.exact) {
        score = MIN_SCORE;
        trim_len = if forced_end == ForceQryTerminus::QryEnd { // as measured from SEQ 5' or 3' end; corresponds to the A-tail base
            i0 + MIN_SCORE as usize - cut_site_offset
        } else {
            SEARCH_SPACE_3 - i0 - cut_site_offset
        };

    // if not found, use smith_waterman to locate the adapter sequence with the RE half-site
    // will have a sequencing error within the 10bp exact search region
    } else {
        let aln = aligner.align(
            &adapter.adapter_re, 
            read_segment, 
            Some(&forced_end), 
            false,
        );
        if aln.score >= MIN_SCORE {
            score = aln.score;
            trim_len = if forced_end == ForceQryTerminus::QryEnd {
                aln.tgt_end0 + 1 - cut_site_offset
            } else {
                SEARCH_SPACE_3 - aln.tgt_start0 - cut_site_offset
            };

        // if not found, use smith_waterman to locate an adapter-gDNA transition without the RE half-site
        // in addition to sequencing errors, this may arise from adapter ligation onto a non-RE gDNA end
        } else if include_re_half_site {
            let aln = aligner.align(
                &adapter.adapter, 
                read_segment, 
                Some(&forced_end), 
                false,
            );
            if aln.score >= MIN_SCORE {
                score = aln.score;
                trim_len = if forced_end == ForceQryTerminus::QryEnd {
                    aln.tgt_end0 + 1
                } else {
                    SEARCH_SPACE_3 - aln.tgt_start0
                };
            }
        }
    }

    // return final results
    (score >= MIN_SCORE, trim_len)

}

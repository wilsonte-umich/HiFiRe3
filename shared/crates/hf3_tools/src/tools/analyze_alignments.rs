//  actions:
//    prepare read alignments for variant analysis, performing steps that don't require RE site data
//      remove unneeded alignment tags
//      calculate the number of reference and read bases in each alignment
//      split unmerged read pairs for downstream processing as single-end reads
//         suffix QNAME with /1 or /2 to make reads unique for independent processing
//      order alignments across each read from 5' to 3'
//      when applicable to targeted libraries:
//         assess whether alignment/read ends matched a target region
//     enforce alignment-level quality filtering
//         used to suppress untrusted junctions arising from aligner issues
//     calculate the traversal delta over all sets of read alignments (across all sets of junctions)
//          used to suppress spurious junctions arising in low quality regions with ~equal inserted and deleted bases
//     add metadata about microhomlogous/inserted bases at SV junctions
//     flag chimeric junctions at SV breakpoint nodes that match:
//         - foldback inversions (due to ONT duplex reads or inversion synthesis)
//         - low quality insertions (ONT only, to suppress follow-on chimeras)
//         - adapter insertions (due to any process that conjoins reads across an adapter)
//      when applicable to RE fragment libraries:
//        extract all sequence endpoints expected to match to potential RE sites later, i.e.:
//          5' endpoints of all sequences, i.e., reads
//          3' endpoints from 3'-adapter-trimmed single-reads or merged read pairs, i.e., molecules sequenced end-to-end
//              *>>>>>>>>>>>o------  (truncated/untrimmed single-read, 3' end not used)
//                *>>>>>>>>>>>>>>>>>>* (merged read pair, 3'-adapter-trimmed ONT or other end-to-end single read)
//            *>>>>>>------<<<<<<* (unmerged paired reads)
//        where "endpoint" means the outermost bases of the sequenced event, ignoring internal supplementary alignments
//            *aln1//aln2//aln3* (aln2 is ignored by this script)
//  expects:
//      source $MODULES_DIR/genome/set_genome_vars.sh
//      source $MODULES_DIR/align/set_alignment_vars.sh
//      source $MODULES_DIR/library/set_library_vars.sh
//  input:
//      SAM stream from reference minimap2 alignment on STDIN
//  output: 
//      modified SAM on STDOUT with custom tags added at end
//          unmerged paired reads are suffixed with /1 or /2 in QNAME for independent processing downstream
//      headerless table as OBSERVED_ENDPOINTS_FILE of observed endpoints with columns:
//        chrom,sitePos1,nObserved
//      where:
//        chrom is the chromosome name per setCanonicalChroms()
//        sitePos1 is the 1-referenced genome position adjusted to infer at least one endpoint's candidate blunt RE site
//                   |*======   where | = RE cleaved bond, + = refPos1, * = sitePos1, side = R
//            ======+|*         and side = L
//                   |*--+===   where clipped (-) but not trimmed bases imply that RE site is offset from refPos1
//            ====+--|*         small clips more likely reflect base differences in contiguous alignments, large clips may represent SV junctions
//        nObserved is the count of all unique endpoints that nominated chrom,sitePos1
//        all sitePos1 are for blunt REs, since read endpoint analysis only applies to ligation libraries

// modules
mod alignment_quality_filters;
mod traversal_delta;
mod chimeric_read_splitter;
mod endpoint_extraction;

// dependencies
use std::error::Error;
use genomex::sam::junction::JunctionType;
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters, COUNTER_SEPARATOR};
use mdi::RecordStreamer;
use genomex::sam::{SamRecord, flag};
use genomex::genome::{Chroms, TargetRegion, TargetRegions};
use crate::formats::hf_tags::*;
use endpoint_extraction::Endpoints;
use alignment_quality_filters::AlnFailure;
use traversal_delta::Traversal;
use chimeric_read_splitter::ChimeraSplitter;
use crate::junctions::JxnFailureFlag;

// Tool structure, for passing data workers to record processing functions.
struct Tool {
    chroms: Chroms,
    targets:       TargetRegions,
    aln_failure:   AlnFailure,
    traversal:     Traversal,
    splitter:      ChimeraSplitter,
    endpoints:     Endpoints,
    incoming_tags: Vec<&'static str>,
    stage_tags:    Vec<&'static str>,
}

// constants
const TOOL: &str = "analyze_alignments";
pub_key_constants!(
    // from environment variables
    SEQUENCING_PLATFORM
    IS_END_TO_END_READ
    READ_PAIR_TYPE
    REJECTING_JUNCTION_RE_SITES
    READ_LENGTH_TYPE
    // MIN_MAPQ
    // CLIP_TOLERANCE
    // HAS_BASE_ACCURACY
    // derived configuration values
    IS_END_TO_END_PLATFORM
    IS_PAIRED_READS
    IS_ONT
    // counter keys
    N_SEQS
    N_READS
    N_ALNS
    N_BLOCKS
    N_BY_GENOME
    N_BY_OUTCOME
    N_READS_OUT
    N_READS_SV
    N_TOTAL_JXNS
    N_JXNS_BY_REASON
    // read outcomes
    USABLE_READ 
    NONCANONICAL_CHROM 
    IS_OFF_TARGET 
    // jxn failure keys
    // JXN_FAIL_NONE // checked by chimera splitter
    JXN_FAIL_TRAVERSAL_DELTA // checked here
    JXN_FAIL_NONCANONICAL
    // JXN_FAIL_FOLDBACK_INV
    // JXN_FAIL_ONT_FOLLOW_ON
    // JXN_FAIL_HAS_ADAPTER
    // JXN_FAIL_SITE_MATCH // not used yet, check later

    // N_JXNS
    // N_BLOCKS
    // N_BAD_TRAVERSAL
);

//     ALN_FAIL_NONE            => 0,
//     ALN_FAIL_MAPQ            => 1,
//     ALN_FAIL_DIVERGENCE      => 2,
//     ALN_FAIL_FLANK_LEN       => 4,
//     ALN_FAIL_AVG_BASE_QUAL   => 8,
//     # -------------
//     JXN_FAIL_NONE            => 0,
//     JXN_FAIL_TRAVERSAL_DELTA => 1,
//     JXN_FAIL_NONCANONICAL    => 2,
//     JXN_FAIL_SITE_MATCH      => 4,
//     JXN_FAIL_FOLDBACK_INV    => 8,
//     JXN_FAIL_ONT_FOLLOW_ON   => 16, 
//     JXN_FAIL_HAS_ADAPTER     => 32,

// tool function called by hf3_tools main()
pub fn stream() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[SEQUENCING_PLATFORM, READ_PAIR_TYPE]);
    cfg.set_bool_env(&[IS_END_TO_END_READ, REJECTING_JUNCTION_RE_SITES]);

    // set derived config values
    cfg.set_bool( IS_END_TO_END_PLATFORM, *cfg.get_bool(IS_END_TO_END_READ));
    cfg.set_bool( IS_PAIRED_READS,        cfg.equals_string(READ_PAIR_TYPE, "paired"));
    cfg.set_bool( IS_ONT,                 cfg.equals_string(SEQUENCING_PLATFORM, "ONT"));

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_SEQS,  "sequences, i.e., single reads or read pairs"),
        (N_READS, "reads"),
        (N_ALNS,  "alignments"),
        (COUNTER_SEPARATOR, ""),
        (N_READS_OUT,  "reads were fully processed"),
        (N_READS_SV,   "fully processed reads with putative structural variants"),
        (N_TOTAL_JXNS, "fully processed reads with putative structural variants"),
        (N_BLOCKS,     "alignment blocks"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_BY_GENOME,      "number of reads by genome"),
        (N_BY_OUTCOME,     "number of reads by outcome"),
        (N_JXNS_BY_REASON, "junction failure counts by reason")
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.print("initializing");

    // build tool support resources
    let mut tool = Tool {
        chroms:        Chroms::new(&mut w.cfg),
        targets:       TargetRegions::new(&mut w),
        aln_failure:   AlnFailure::new(&mut w),
        traversal:     Traversal::new(&mut w),
        splitter:      ChimeraSplitter::new(&mut w),
        endpoints:     Endpoints::new(&mut w),
        incoming_tags: StageTags::Alignment.tags_after_stage(),
        stage_tags:    StageTags::AlignmentAnalysis.tag_added_by_stage(),
    };

    // process groups of SAM records for each read in a stream
    w.log.print("processing reads from SAM records");
    let mut rs = RecordStreamer::new();
    rs
        .comment(b'@')
        .no_trim()
        .flexible()
        .group_by_in_place_serial(|alns: &mut Vec<SamRecord>| record_parser(
            alns, 
            &mut w,
            &mut tool,
        ), &["qname"]);

    // report counter values
    Endpoints::write(&mut w, &mut tool)?;
    w.ctrs.print_all();
    Ok(())
}

// process all alignments for one (paired) read
fn record_parser(
    alns: &mut Vec<SamRecord>, 
    w: &mut Workflow,
    tool: &mut Tool,
) -> Result<Vec<usize>, Box<dyn Error>> {
    w.ctrs.increment(N_SEQS);

    // remove unneeded alignment tags
    // calculate the number of reference and read bases in each alignment
    // split alignment by read number and assess pairing status
    let mut read1_alns: Vec<&mut SamRecord> = Vec::new();
    let mut read2_alns: Vec<&mut SamRecord> = Vec::new();
    for aln in alns.iter_mut() {
        aln.tags.retain(&tool.incoming_tags);
        aln.tags.tags.push(format!("{}{}", N_READ_BASES, aln.get_aligned_size()));
        aln.tags.tags.push(format!("{}{}", N_REF_BASES,  aln.get_ref_span()));
        if aln.check_flag_all(flag::IS_PAIRED + flag::SECOND_IN_PAIR) {
            read2_alns.push(aln);
        } else {
            read1_alns.push(aln);
        }
    }
    let is_unmerged_pair: bool = !read2_alns.is_empty();

    // sort alignments by query start position
    // must do upstream of paired_node5 assignment
    read1_alns.sort_by_key(|aln| aln.get_query_start0());
    if is_unmerged_pair {
        read2_alns.sort_by_key(|aln| aln.get_query_start0());
    }

    // process each read's alignments
    let paired_node5 = if is_unmerged_pair { &get_paired_node5(tool, read2_alns[0]) } else { "" };
    process_read(1, &mut read1_alns, w, tool, is_unmerged_pair, paired_node5)?;
    if is_unmerged_pair {
        let paired_node5 = &get_paired_node5(tool, read1_alns[0]);
        process_read(2, &mut read2_alns, w, tool, is_unmerged_pair, paired_node5)?;
    }

    // SAM record processed successfully ready for printing to STDOUT
    Ok((0..alns.len()).collect())
}

// get the 5' alignment node of the paired read
fn get_paired_node5(tool: &Tool, aln: &SamRecord) -> String {
    if aln.check_flag_any(flag::UNMAPPED) || !tool.chroms.is_canonical(&aln.rname) { 
        return "*".to_string();
    }
    let (pos1, strand) = if aln.check_flag_any(flag::REVERSE) { 
        (aln.get_end1(), '-') 
    } else { 
        (aln.pos1, '+') 
    };
    format!("{}:{}{}", aln.rname, pos1, strand)
}

// process the alignments for one read (of a pair)
fn process_read(
    read_n: usize, 
    alns: &mut Vec<&mut SamRecord>, 
    w: &mut Workflow,
    tool: &mut Tool,
    is_unmerged_pair: bool,
    paired_node5: &str,
) -> Result<(), Box<dyn Error>> {
    w.ctrs.increment(N_READS);
    w.ctrs.add_to(N_ALNS, alns.len());

    // adjust QNAMEs to indicate unmerged read number, i.e., to be unique per paired read
    // from here forward, unmerged paired reads are treated as fixed-length single reads
    // junctions in anomalous gaps are unverifiable and not considered for rare SV analysis
    // pass the paired 5' read end for proper outer node assignment
    if is_unmerged_pair {
        for aln in alns.iter_mut() {
            aln.qname = format!("{}/{}", aln.qname, read_n);
            aln.tags.tags.push(format!("{}{}", UNMERGED_NODE5, paired_node5));
        }
    }

    // stop processing reads further whose 5' alignment is not on a canonical chromosome
    // they are either concretely or effectively unmapped
    if alns[0].check_flag_any(flag::UNMAPPED) || !tool.chroms.is_canonical(&alns[0].rname) { 
        alns[0].seq       = "*".to_string();
        alns[0].qual.qual = "*".to_string();
        w.ctrs.increment_keyed(N_BY_OUTCOME, NONCANONICAL_CHROM);
        return Ok(()); 
    }

    // count reads by genome to report on mixed species samples
    w.ctrs.increment_keyed(N_BY_GENOME, &tool.chroms.get_genome_suffix(&alns[0].rname));

    // determine the target status of the read, which is always "on target" for non-targeted libraries
    let is_on_target = tool.targets.set_aln_target_classes(
        alns, w, TARGET_CLASS, IS_ON_TARGET
    );

    // stop processing off-target reads, we have all required counts and they won't be used for variant calling
    if !is_on_target {
        alns[0].seq       = "*".to_string();
        alns[0].qual.qual = "*".to_string();
        w.ctrs.increment_keyed(N_BY_OUTCOME, IS_OFF_TARGET);
        return Ok(());
    }

    // add alignment-level metadata
    // assess alignment-level quality filtering (not a jxn rejection yet, but used to filter SVs downstream)
    let n_alns = alns.len();
    let read_has_jxn = n_alns > 1;
    tool.aln_failure.set_aln_failure_flag(w, alns[0], read_has_jxn);
    let mut block_n = 1;
    alns[0].tags.tags.push(format!("{}{}", BLOCK_N, block_n));
    if read_has_jxn {
        w.ctrs.increment(N_READS_SV);
        w.ctrs.add_to(N_TOTAL_JXNS, n_alns - 1);
        
        // assess traversal delta for every pair of alignments flanking every possible read sub-path
        // fail all junctions within a failed path
        // once a junction fails, it stays failed even if it passes a different pair of alignments
        let mut failed_traversal_delta = vec![false; n_alns];
        for aln2_i in 1..n_alns {
            for aln1_i in 0..aln2_i {
                if tool.traversal.failed_delta(alns[aln1_i], alns[aln2_i]) {
                    for fail_i in (aln1_i + 1)..=aln2_i {
                        failed_traversal_delta[fail_i] = true;
                    }
                }
            }
        }

        // set block_n based on failed_traversal_delta
        // check junctions for traversal delta failure and other criteria
        for aln2_i in 1..n_alns {
            tool.aln_failure.set_aln_failure_flag(w, alns[aln2_i], read_has_jxn);

            let jxn = SamRecord::get_junction_metadata(alns[aln2_i - 1], alns[aln2_i]);

            if failed_traversal_delta[aln2_i] {
                w.ctrs.increment_keyed(N_JXNS_BY_REASON, JXN_FAIL_TRAVERSAL_DELTA);
                alns[aln2_i].tags.tags.push(format!("{}{}", BLOCK_N, block_n));
                alns[aln2_i - 1].tags.tags.push(format!("{}{}", JXN_FAILURE_FLAG_TMP, JxnFailureFlag::TraversalDelta as u8));
            } else {
                block_n += 1; // increment alignment blockN on every non-failed junction
                alns[aln2_i].tags.tags.push(format!("{}{}", BLOCK_N, block_n));
                let jxn_failure_flag = if 
                    !tool.chroms.is_canonical(&alns[aln2_i - 1].rname) || 
                    !tool.chroms.is_canonical(&alns[aln2_i].rname) 
                {
                    // check non-canonical chromosome junctions here to avoid mutable borrow issues
                    w.ctrs.increment_keyed(N_JXNS_BY_REASON, JXN_FAIL_NONCANONICAL);
                    JxnFailureFlag::Noncanonical
                } else {
                    // check other junction failure reasons (get_jxn_failure_flag handles N_JXNS_BY_REASON counter increments)
                    tool.splitter.get_jxn_failure_flag(w, alns[aln2_i - 1], alns[aln2_i], &jxn)
                };
                alns[aln2_i - 1].tags.tags.push(format!("{}{}", JXN_FAILURE_FLAG_TMP, jxn_failure_flag as u8));
            }
            alns[aln2_i - 1].tags.tags.push(format!("{}{}", JXN_TYPE,   jxn.jxn_type as u8));
            alns[aln2_i - 1].tags.tags.push(format!("{}{}", ALN_OFFSET, jxn.alignment_offset));
            alns[aln2_i - 1].tags.tags.push(format!("{}{}", JXN_BASES,  jxn.jxn_bases));
        }
    
    // drop SEQ and QUAL on invariant junctions to save disk space
    // TODO: retain these if rare SNV analysis is intended downstream?
    } else {
        alns[0].seq       = "*".to_string();
        alns[0].qual.qual = "*".to_string();
    }

    // fill in null junction metadata on the 3' most aln, from JXN_FAILURE_FLAG_TMP to JXN_BASES
    let i_3 = n_alns - 1;
    alns[i_3].tags.tags.push(format!("{}{}", JXN_FAILURE_FLAG_TMP, JxnFailureFlag::None as u8));
    alns[i_3].tags.tags.push(format!("{}{}", JXN_TYPE,   JunctionType::Proper as u8));
    alns[i_3].tags.tags.push(format!("{}{}", ALN_OFFSET, 0));
    alns[i_3].tags.tags.push(format!("{}{}", JXN_BASES,  "*"));

    // extract endpoints as possible RE sites
    Endpoints::extract(read_n, alns, w, tool, is_unmerged_pair)?;

    // increment usable read counter and return
    w.ctrs.increment(N_READS_OUT);
    w.ctrs.increment_keyed(N_BY_OUTCOME, USABLE_READ);
    Ok(())
}

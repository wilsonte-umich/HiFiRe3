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
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use mdi::RecordStreamer;
use genomex::sam::{SamRecord, flag};
use genomex::genome::{Chroms, TargetRegions};
use genomex::sam::nullable::*;
use crate::formats::hf_tags::*;
use crate::junctions::JxnFailureFlag;
use alignment_quality_filters::{AlnFailure, N_ALNS_BY_REASON, AVG_BASE_QUAL};
use traversal_delta::Traversal;
use chimeric_read_splitter::{ChimeraSplitter, ADAPTER_SCORES, N_CHIMERIC};
use endpoint_extraction::{Endpoints, N_ENDS, N_HIGH_QUAL, N_ENDPOINTS};

// Tool structure, for passing data workers to record processing functions.
struct Tool {
    chroms: Chroms,
    targets:       TargetRegions,
    aln_failure:   AlnFailure,
    traversal:     Traversal,
    splitter:      ChimeraSplitter,
    endpoints:     Endpoints,
    incoming_tags: Vec<&'static str>,
    has_base_accuracy: bool,
}

// constants
const TOOL: &str = "analyze_alignments";
pub_key_constants!(
    // from environment variables
    SEQUENCING_PLATFORM
    IS_END_TO_END_READ
    READ_PAIR_TYPE
    HAS_BASE_ACCURACY
    // derived configuration values
    IS_END_TO_END_PLATFORM
    IS_PAIRED_READS
    IS_ONT
    // counter keys
    N_SEQS
    N_READS
    N_ALNS
    N_BLOCKS
    N_READS_OUT
    N_READS_SV
    N_TOTAL_JXNS
    N_BY_GENOME
    N_READS_BY_OUTCOME
    N_JXNS_BY_REASON
    // read outcomes
    UNMAPPED_READ
    NONCANONICAL_CHROM 
    IS_OFF_TARGET 
    USABLE_READ 
    // jxn failure keys
    JXN_FAIL_NONE
    JXN_FAIL_TRAVERSAL_DELTA
    JXN_FAIL_NONCANONICAL
    JXN_FAIL_FOLDBACK_INV
    JXN_FAIL_ONT_FOLLOW_ON
    JXN_FAIL_HAS_ADAPTER
    JXN_FAIL_SITE_MATCH // not used yet, check later
    JXN_FAIL_STEM_LENGTH
    JXN_FAIL_UPSTREAM
);

// tool function called by hf3_tools main()
pub fn stream() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[SEQUENCING_PLATFORM, READ_PAIR_TYPE]);
    cfg.set_bool_env(&[IS_END_TO_END_READ, HAS_BASE_ACCURACY]);

    // set derived config values
    cfg.set_bool( IS_END_TO_END_PLATFORM, *cfg.get_bool(IS_END_TO_END_READ));
    cfg.set_bool( IS_PAIRED_READS,        cfg.equals_string(READ_PAIR_TYPE, "paired"));
    cfg.set_bool( IS_ONT,                 cfg.equals_string(SEQUENCING_PLATFORM, "ONT"));

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_SEQS,       "input sequences, i.e., single reads or read pairs"),
        (N_READS,      "input reads"),
        (N_ALNS,       "input alignments"),
        (N_READS_OUT,  "fully processed on-target reads"),
        (N_READS_SV,   "fully processed reads with putative structural variants"),
        (N_TOTAL_JXNS, "junction in fully processed reads"),
        (N_BLOCKS,     "alignment blocks"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_READS_BY_OUTCOME, "number of reads by outcome"),
        (N_BY_GENOME,        "number of reads by genome (canonical chroms only)"),
        (N_JXNS_BY_REASON,   "junction failure counts by reason")
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
        has_base_accuracy: *w.cfg.get_bool(HAS_BASE_ACCURACY),
    };

    // initialize the record streamer
    let mut rs = RecordStreamer::new();
    rs
        .comment(b'@') // pass SAM header lines
        .no_trim()
        .flexible();   // to support variable numbers of tags

    // process groups of SAM records for each read in a stream
    // use different record_parsers for single vs. paired-end reads
    w.log.print("processing reads from SAM records");
    if *w.cfg.get_bool(IS_PAIRED_READS) {
        rs.group_by_in_place_serial(
            |alns: &mut [SamRecord]| parse_paired_end(alns, &mut w, &mut tool), 
            &["qname"]
        );
    } else {
        rs.group_by_in_place_serial(
            |alns: &mut [SamRecord]| parse_single_read(alns, &mut w, &mut tool), 
            &["qname"]
        );
    }; 

    // write the aggregated endpoint counts
    tool.endpoints.write(&mut w, &tool.chroms)?;

    // report counter values
    w.ctrs.print_grouped(&[
        &[N_SEQS, N_READS, N_ALNS],
        &[N_READS_BY_OUTCOME],
        &[N_BY_GENOME],
        &[N_READS_OUT, N_READS_SV, N_TOTAL_JXNS, N_BLOCKS, N_CHIMERIC],
        &[N_ALNS_BY_REASON],
        &[AVG_BASE_QUAL],
        &[N_JXNS_BY_REASON],
        &[ADAPTER_SCORES],
        &[N_ENDS, N_HIGH_QUAL, N_ENDPOINTS],
    ]);
    Ok(())
}

// process all alignments for a paired-end read pair
fn parse_paired_end(
    alns: &mut [SamRecord], 
    w: &mut Workflow,
    tool: &mut Tool,
) -> Result<Vec<usize>, Box<dyn Error>> {

    // short-circuit merged read pairs by treating them as single reads
    // regardless of how many alignments the merged read has
    if alns[0].get_tag_value(FASTP_MERGE).is_some(){
        return parse_single_read(alns, w, tool);
    }
    w.ctrs.increment(N_SEQS);

    // adjust QNAMEs to indicate unmerged read number, i.e., to be unique per paired read
    // from here forward, unmerged paired reads are treated as fixed-length single reads
    // junctions in anomalous gaps are unverifiable and not considered for rare SV analysis
    let mut n_read1: usize  = 0;
    for aln in alns.iter_mut() {
        let read_n: usize = if aln.check_flag_all(flag::IS_PAIRED + flag::SECOND_IN_PAIR) { 2 } else { 1 };
        if read_n == 1 { n_read1 += 1; }
        aln.qname = format!("{}/{}", aln.qname, read_n);
    }

    // sort unmerged paired reads by read number then by query start position
    alns.sort_by_key(|aln| {(
        if aln.qname.ends_with("/2") { 2 } else { 1 }, 
        aln.get_query_start0()
    )});
    
    // process each read's 5'-3' sorted alignments
    let paired_outer_node = alns[n_read1].pack_signed_node_aln(5, &tool.chroms);
    process_read(1, &mut alns[0..n_read1], w, tool, true, paired_outer_node)?;
    let paired_outer_node = alns[0].pack_signed_node_aln(5, &tool.chroms);
    process_read(2, &mut alns[n_read1..], w, tool, true, paired_outer_node)?;
    
    // print all alignments for both paired reads, regardless of outcome
    Ok((0..alns.len()).collect())
}

// process all alignments for a single read or merged read pair
fn parse_single_read(
    alns: &mut [SamRecord], 
    w: &mut Workflow,
    tool: &mut Tool,
) -> Result<Vec<usize>, Box<dyn Error>> {
    w.ctrs.increment(N_SEQS);

    // sort multiple alignments by query start position
    let n_alns = alns.len();
    if n_alns > 1 {
        alns.sort_by_key(|aln| aln.get_query_start0());
    }

    // process the 5'-3' sorted alignments
    process_read(1, alns, w, tool, false, 0)?;

    // print all alignments for the read, regardless of outcome
    Ok((0..n_alns).collect())
}

// process the alignments for one read (of a pair)
fn process_read(
    read_n: usize, 
    alns: &mut [SamRecord], 
    w: &mut Workflow,
    tool: &mut Tool,
    is_unmerged_pair: bool,
    paired_outer_node: isize,
) -> Result<(), Box<dyn Error>> {
    w.ctrs.increment(N_READS);
    let n_alns = alns.len();
    w.ctrs.add_to(N_ALNS, alns.len());

    // remove unneeded alignment tags, mate information
    // pass the paired 5' read end for proper outer node assignment
    for aln in alns.iter_mut() {
        aln.set_to_null(&[RNEXT, PNEXT, TLEN]);
        aln.tags.retain(&tool.incoming_tags);
        if is_unmerged_pair { aln.tags.tags.push(format!("{}{}", PAIRED_OUTER_NODE, paired_outer_node)); }
    }

    // stop processing unmapped reads
    if alns[0].check_flag_any(flag::UNMAPPED) { 
        alns[0].set_tag_value(INSERT_SIZE, &alns[0].seq.len()); // set insert size for unmapped reads
        alns[0].set_to_null(&[SEQ, QUAL]); // then drop the sequence and quality
        w.ctrs.increment_keyed(N_READS_BY_OUTCOME, UNMAPPED_READ);
        return Ok(()); 
    }

    // stop processing reads further when the 5' alignment is not on a canonical chromosome
    // they are either concretely or effectively unmapped
    // regard them as off-target, so set tag READ_IS_OFF_TARGET
    if !tool.chroms.is_canonical(&alns[0].rname) { 
        for aln in alns.iter_mut() { 
            aln.set_to_null(&[SEQ, QUAL]); 
            aln.tags.drop(&[DIFFERENCE_STRING]);
            aln.set_tag_value(READ_IS_OFF_TARGET, &1);
        }
        w.ctrs.increment_keyed(N_READS_BY_OUTCOME, NONCANONICAL_CHROM);
        return Ok(()); 
    }

    // count reads by genome to report on mixed species samples
    w.ctrs.increment_keyed(N_BY_GENOME, &tool.chroms.get_genome_suffix(&alns[0].rname));

    // determine the target status of the read, which is always "on target" for non-targeted libraries
    let is_on_target = tool.targets.set_aln_target_classes(
        alns, w, TARGET_CLASS, READ_IS_OFF_TARGET
    );

    // stop processing off-target reads, we have all required counts and they won't be used for variant calling
    if !is_on_target {
        for aln in alns.iter_mut() { 
            aln.set_to_null(&[SEQ, QUAL]); 
            aln.tags.drop(&[DIFFERENCE_STRING]) // READ_IS_OFF_TARGET was set by set_aln_target_classes
        }
        w.ctrs.increment_keyed(N_READS_BY_OUTCOME, IS_OFF_TARGET);
        return Ok(());
    }

    // assess alignment-level quality filtering on the first alignment
    // not a jxn rejection yet, but used to filter SVs downstream
    let aln5_i = 0_usize;
    let read_has_jxn = n_alns > 1;
    let mut block_n = 1;
    tool.aln_failure.set_aln_failure_flag(w, &mut alns[aln5_i], read_has_jxn);

    // continue assessing alignment and junction-level quality filtering
    // for all alignements beyond the first
    // record junction information on the alignment 5'-proximal to the junction
    if read_has_jxn {
        w.ctrs.increment(N_READS_SV);
        w.ctrs.add_to(N_TOTAL_JXNS, n_alns - 1);
        
        // assess traversal delta for every pair of alignments flanking every possible read sub-path
        // fail all junctions within a failed path
        // once a junction fails, it stays failed even if it passes a different pair of alignments
        let mut failed_traversal_delta = vec![false; n_alns];
        for aln3_i in 1..n_alns {
            for aln5_i in 0..aln3_i {
                if tool.traversal.failed_delta(&alns[aln5_i], &alns[aln3_i]) {
                    for fail_i in (aln5_i + 1)..=aln3_i {
                        failed_traversal_delta[fail_i] = true;
                    }
                }
            }
        }

        // set block_n based on failed_traversal_delta
        // check junctions for traversal delta failure and other criteria
        for aln3_i in 1..n_alns {
            let aln5_i = aln3_i - 1;
            tool.aln_failure.set_aln_failure_flag(w, &mut alns[aln3_i], read_has_jxn);
            let jxn = SamRecord::get_junction(&alns[aln5_i], &alns[aln3_i], &tool.chroms);
            let jxn_failure_flag = if failed_traversal_delta[aln3_i] {
                w.ctrs.increment_keyed(N_JXNS_BY_REASON, JXN_FAIL_TRAVERSAL_DELTA);
                alns[aln3_i].tags.tags.push(format!("{}{}", BLOCK_N, block_n));
                JxnFailureFlag::TraversalDelta
            } else {
                block_n += 1; // increment alignment blockN on every non-failed junction
                alns[aln3_i].tags.tags.push(format!("{}{}", BLOCK_N, block_n));
                if 
                    !tool.chroms.is_canonical(&alns[aln5_i].rname) || 
                    !tool.chroms.is_canonical(&alns[aln3_i].rname) 
                {
                    // check non-canonical chromosome junctions here to avoid mutable borrow issues
                    w.ctrs.increment_keyed(N_JXNS_BY_REASON, JXN_FAIL_NONCANONICAL);
                    JxnFailureFlag::Noncanonical
                } else {
                    // check other junction failure reasons (get_jxn_failure_flag handles N_JXNS_BY_REASON counter increments)
                    tool.splitter.get_jxn_failure_flag(w, alns, aln5_i, aln3_i, &jxn) 
                }
            };
            alns[aln5_i].tags.tags.push(format!("{}{}", JUNCTION, jxn.serialize()));
            if jxn_failure_flag != JxnFailureFlag::None {
                alns[aln5_i].tags.tags.push(format!("{}{}", JXN_FAILURE_FLAG_INIT, jxn_failure_flag as u8));
            }
        }
    
        // now that junctions are described, drop SEQ and QUAL on all alignments except the 5' alignment
        for aln3_i in 1..n_alns {
            alns[aln3_i].set_to_null(&[SEQ, QUAL]);
        }

    // to save disk space
    //  - drop SEQ and QUAL on reads with no junctions
    //  - drop cs tag unless platform HAS_BASE_ACCURACY, so has short tags that may be used for SNV calling
    } else {
        alns[aln5_i].set_to_null(&[SEQ, QUAL]);
        if !tool.has_base_accuracy {
            alns[aln5_i].tags.drop(&[DIFFERENCE_STRING])
        }
    }

    // extract endpoints as possible RE sites
    tool.endpoints.extract(read_n, alns, w, &tool.chroms, is_unmerged_pair)?;

    // increment usable read counter
    w.ctrs.add_to(N_BLOCKS, block_n);
    w.ctrs.increment(N_READS_OUT);
    w.ctrs.increment_keyed(N_READS_BY_OUTCOME, USABLE_READ);

    // done
    Ok(())
}

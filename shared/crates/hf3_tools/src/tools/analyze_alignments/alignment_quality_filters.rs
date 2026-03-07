//! action:
//!   check whether alignments pass alignment-level quality filters:
//!       all alignments, even single, non-SV alignments per read:
//!           minimum MAPQ
//!           maximium gap-corrected base divergence
//!       alignment in SV reads only:
//!           minimum alignment length, measured as reference base span
//!           (ONT only) minimum average base quality across alignment
//!               with Dorado, individual base qualities range from Q0 to Q50
//!               see script output for the range of average base qualities observed over SV flanking alignments
//!   return a bit-encoded alignment failure flag to suppress untrusted SV junctions that would include a failed alignment

// dependencies
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use genomex::sam::SamRecord;
use crate::formats::hf3_tags::ALN_FAILURE_FLAG;

// constants
pub_key_constants!{
    // from environment variables
    MIN_MAPQ
    MAX_DIVERGENCE
    MIN_FLANK_LEN
    MIN_AVG_BASE_QUAL
    // counters
    N_ALNS_BY_REASON
    AVG_BASE_QUAL
    // alignment failure keys
    ALN_FAIL_NONE
    ALN_FAIL_MAPQ
    ALN_FAIL_DIVERGENCE
    ALN_FAIL_FLANK_LEN
    ALN_FAIL_AVG_BASE_QUAL
}
const AVG_BASE_QUAL_BIN_SIZE: usize = 5;

// alignment failure flags
#[repr(usize)]
pub enum AlnFailureFlag {
    // None        = 0, // implied by absence of alignment failure tag
    Mapq        = 1,
    Divergence  = 2,
    FlankLen    = 4,
    AvgBaseQual = 8,
}

/// Traversal structure for assessing low quality internal spans
/// with likely erroneous alignments.
pub struct AlnFailure {
    min_mapq:          u8,
    max_divergence:    f64,
    min_flank_len:     u32,
    min_avg_base_qual: f64,
}
impl AlnFailure {

    /// Initialize an alignment failure checker.
    pub fn new(w: &mut Workflow) -> AlnFailure {
        w.cfg.set_u8_env(&[MIN_MAPQ]);
        w.cfg.set_f64_env(&[MAX_DIVERGENCE, MIN_AVG_BASE_QUAL]);
        w.cfg.set_u32_env(&[MIN_FLANK_LEN]);
        w.ctrs.add_keyed_counters(&[
            (N_ALNS_BY_REASON, "alignments failure counts by reason")
        ]);
        w.ctrs.add_indexed_counters(&[
            (AVG_BASE_QUAL, "qual_bin", 10, "SV alignment average base quality bin counts (multiply by 5 to get QUAL)"),
        ]);
        AlnFailure{
            min_mapq:          *w.cfg.get_u8(MIN_MAPQ),
            max_divergence:    *w.cfg.get_f64(MAX_DIVERGENCE),
            min_flank_len:     *w.cfg.get_u32(MIN_FLANK_LEN),
            min_avg_base_qual: *w.cfg.get_f64(MIN_AVG_BASE_QUAL),
        }
    }

    // Check the quality of an individual alignment, returning an alignment failure flag.
    pub fn set_aln_failure_flag(&self, w: &mut Workflow, aln: &mut SamRecord, read_has_jxn: bool) -> bool {

        // rejection criteria are enforced sequentially in order of efficiency
        // i.e., frequent, easy rejections are checked first
        // later criteria are not check if an earlier criterion already failed

        // criteria enforced on all alignments, even single, non-SV alignments
        if aln.mapq < self.min_mapq {
            w.ctrs.increment_keyed(N_ALNS_BY_REASON, ALN_FAIL_MAPQ);
            aln.tags.tags.push(format!("{}{}", ALN_FAILURE_FLAG, AlnFailureFlag::Mapq as u8));
            return false;
        }
        let de = aln.get_tag_value_parsed("de").unwrap_or(0.0);
        if de > self.max_divergence {
            w.ctrs.increment_keyed(N_ALNS_BY_REASON, ALN_FAIL_DIVERGENCE);
            aln.tags.tags.push(format!("{}{}", ALN_FAILURE_FLAG, AlnFailureFlag::Divergence as u8));
            return false;
        }

        // criteria only enforced when reads have SV junctions, i.e., multiple alignments
        if read_has_jxn {
            if aln.get_ref_span() < self.min_flank_len {
                w.ctrs.increment_keyed(N_ALNS_BY_REASON, ALN_FAIL_FLANK_LEN);
                aln.tags.tags.push(format!("{}{}", ALN_FAILURE_FLAG, AlnFailureFlag::FlankLen as u8));
                return false;
            }
            if self.min_avg_base_qual > 0.0 { // this slow check only performed on platforms that need it (ONT, PacBio)
                let avg_base_qual = aln.get_avg_qual_aln();
                let base_qual_bin_i = (avg_base_qual / AVG_BASE_QUAL_BIN_SIZE as f64).round() as usize;
                w.ctrs.increment_indexed(AVG_BASE_QUAL, base_qual_bin_i);
                if avg_base_qual < self.min_avg_base_qual {
                    w.ctrs.increment_keyed(N_ALNS_BY_REASON, ALN_FAIL_AVG_BASE_QUAL);
                    aln.tags.tags.push(format!("{}{}", ALN_FAILURE_FLAG, AlnFailureFlag::AvgBaseQual as u8));
                    return false;
                }
            }
        }

        // alignment passed all criteria
        w.ctrs.increment_keyed(N_ALNS_BY_REASON, ALN_FAIL_NONE);
        true
    }

}

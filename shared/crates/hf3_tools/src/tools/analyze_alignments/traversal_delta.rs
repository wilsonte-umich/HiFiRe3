// action:
//   calculate the reference vs. read traversal delta between all sets of alignments (across all sets of junctions)
//   use calculated deltas and MIN_TRAVERSAL_DELTA to set block number of all alignments
//       any set of alignments where traversalDelta < MIN_TRAVERSAL_DELTA are assigned the same block number
//   junctions between alignments with the same block number may be considered artifactual during SV analysis
//   block numbering is done before alignment quality truncation when all read alignments are still present
//       unlike other SV metrics, block numbering can depend on alignments distant from the two flanking a junction if nAln > 2
//       we want to ensure the even lower quality alignments have an opportunity to question a read's traversal delta
//   traversal deltas will fail and suppress false junctions when:
//       a large, low quality base string is called as a deletion SV with an equally large insertion
//       a similar base string is falsely aligned to another genomic location when it should have aligned inline
//           =========*==========*========= reference locus 1
//           ||||||||||          ||||||||||
//           ---------*~~~~~~~~~~*--------- query read, where lower quality ~~ bases are either an unaligned insertion or improperly aligned 
//                     ||||||||||
//           ============================== reference locus 2
//       notice that the traversal on reference and query are ~the same from * to *
//       contrast the above with a true deletion with some smaller number of bases inserted at the junction
//           =========*==========*========= reference locus 1
//           ||||||||||          ||||||||||
//           ---------*    ~~    *--------- query read, where ~~ bases are a junction insertion derived from the joining mechanism
//    another way of stating it is that 1-4 is co-linear on query and reference if intervening 1-2 and 3-4 junctions are artifactual
//    such that the first and last segments would not call a SV if the read were aligned end-to-end to reference
//           =========1==========4=========
//           ||||||||||          ||||||||||
//           ---------1~~~~~~~~~~4---------
//                     ||||||||||
//           ==========2========3=========

// dependencies
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use genomex::sam::{SamRecord, flag};

// constants
pub_key_constants!{
    // from environment variables
    MIN_TRAVERSAL_DELTA
}

/// Traversal structure for assessing low quality internal spans
/// with likely erroneous alignments.
pub struct Traversal {
    min_traversal_delta: u32,
}
impl Traversal {

    /// Initialize traversal delta checks.
    pub fn new(w: &mut Workflow) -> Traversal {
        w.cfg.set_u32_env(&[MIN_TRAVERSAL_DELTA]);
        Traversal {
            min_traversal_delta: *w.cfg.get_u32(MIN_TRAVERSAL_DELTA),
        }
    }

    /// Report `false` if traversal delta is consistent with a true SV, `true` if traversal delta fails.
    /// ``````
    /// =========1==========2========= reference locus flanking the query span
    /// ||||||||||          ||||||||||
    /// ---------1~~~~~~~~~~2--------- query read, where low quality ~~ bases are either
    ///           ||||||||||              an unmapped insertion or improperly aligned
    /// ==========2========1========== a chain of 0 to N false alignments
    /// ```
    pub fn failed_delta(&self, aln5: &SamRecord, aln3: &SamRecord) -> bool {
        if aln5.rname != aln3.rname { return false; } // translocation, always a passing traversal delta
        let is_reverse5 = aln5.check_flag_any(flag::REVERSE);
        let is_reverse3 = aln3.check_flag_any(flag::REVERSE);
        if is_reverse5 != is_reverse3 { return false; } // inversion, always a passing traversal delta
        let ref_pos1_aln5_end3 = if is_reverse5 { aln5.pos1 } else { aln5.get_end1() }; // i.e., the inner node positions flanking the junction(s)
        let ref_pos1_aln3_end5 = if is_reverse3 { aln3.get_end1() } else { aln3.pos1 }; 
        let qry_pos1_aln5_end3 = aln5.get_query_end1();
        let qry_pos1_aln3_end5 = aln3.get_query_start0() + 1;
        let ref_traversal = if is_reverse5 { ref_pos1_aln5_end3 - ref_pos1_aln3_end5 } 
                                             else { ref_pos1_aln3_end5 - ref_pos1_aln5_end3 };
        let qry_traversal = qry_pos1_aln3_end5 - qry_pos1_aln5_end3;
        (ref_traversal).abs_diff(qry_traversal) < self.min_traversal_delta
    }
}

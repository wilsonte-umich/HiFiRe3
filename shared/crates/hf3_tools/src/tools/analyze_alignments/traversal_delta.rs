// action:
//   calculate the reference vs. read traversal delta between all sets of alignments (across all sets of junctions)
//   use calculated deltas and MIN_TRAVERSAL_DELTA to set BLOCK_N of all alignments
//       any set of alignments where traversalDelta < MIN_TRAVERSAL_DELTA are assigned the same BLOCK_N
//   junctions between alignments with the same BLOCK_N may be considered artifactual during SV analysis
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
    delta: usize,
}
impl Traversal {

    /// Initialize traversal delta checks.
    pub fn new(w: &mut Workflow) -> Traversal {
        w.cfg.set_usize_env(&[MIN_TRAVERSAL_DELTA]);
        Traversal {
            delta: *w.cfg.get_usize(MIN_TRAVERSAL_DELTA),
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
    pub fn failed_delta(&self, aln1: &SamRecord, aln2: &SamRecord) -> bool {
        if aln1.rname != aln2.rname { return false; } // translocation, always a passing traversal delta
        let is_reverse1 = aln1.check_flag_any(flag::REVERSE);
        let is_reverse2 = aln2.check_flag_any(flag::REVERSE);
        if is_reverse1 != is_reverse2 { return false; } // inversion, always a passing traversal delta
        let ref_pos1_1 = if is_reverse1 { aln1.pos1 } else { aln1.get_end1() }; // i.e., the inner node positions flanking the junction(s)
        let ref_pos1_2 = if is_reverse2 { aln2.get_end1() } else { aln2.pos1 }; 
        let qry_pos1_1 = aln1.get_query_end1();
        let qry_pos1_2 = aln2.get_query_start0() + 1;
        let ref_traversal = if is_reverse1 { ref_pos1_1 - ref_pos1_2 } else { ref_pos1_2 - ref_pos1_1 };
        let qry_traversal = qry_pos1_2 - qry_pos1_1;
        (ref_traversal).abs_diff(qry_traversal) < self.delta
    }
}

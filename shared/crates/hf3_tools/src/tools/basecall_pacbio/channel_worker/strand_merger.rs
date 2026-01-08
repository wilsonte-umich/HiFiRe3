//! Support for comparing and mergin two PacBio strands into a single 
//! duplex consensus read.

// dependencies
use crossbeam::channel::Sender;
use genomex::sequence::Alignment;
use crate::tools::basecall_pacbio::channel_kinetics::kinetics;
use super::{StrandPair, KineticsInstance};

// constants
const PERFECT_MATCH: u8    = 0;
const HAS_END_CLIP: u8     = 1;
const HAS_INDEL: u8        = 2;
const HAS_SUBSTITUTION: u8 = 4;

/// StrandMerger structure for merging two PacBio strands.
pub struct StrandMerger{
    ident_len:  usize,
    del_buffer: String,
    outcomes:   u8,
}
impl StrandMerger {
    /// Create a new StrandMerger structure.
    pub fn new() -> Self {
        StrandMerger { 
            ident_len:  0,
            del_buffer: String::with_capacity(100),
            outcomes:   PERFECT_MATCH,
        }
    }
    /// Merge two PacBio strands into a single consensus read, masking
    /// positions where the strands disagree to N.
    /// 
    /// Create a difference tag similar to minimap2 cs tag to record the
    /// strand differences and the kinetic parameters at divergent positions.
    /// 
    /// Between-strand indels are always recorded as insertions to denote
    /// the possible presence of the base (otherwise it would disappear from
    /// the output seq).
    pub fn merge_strands(
        &mut self,
        strand_pair: &StrandPair,
        aln: &Alignment, // qry in this aligment is the reverse-complemented prev_strand from strand_pair
        tx_kinetics: &Sender<KineticsInstance>,
    ) -> (String, Vec<u8>, String, Vec<u16>, u8) {

        // initialize merger state
        let mut qry_offset = aln.qry_start0;
        let n_qry_bases = strand_pair.prev.seq.len();
        let n_tgt_bases = strand_pair.this.seq.len();
        let mut seq    = String::with_capacity(n_tgt_bases * 2); // the new merged sequence, with Ns at disagreement positions
        let mut dd_tag = String::with_capacity(n_tgt_bases / 4); // the strand differences tag, recording bases at N positions in seq
        let mut sk_tag: Vec<u16> = vec![];
        self.ident_len  = 0;     // length of a current identical run; committed to dd_tag when a difference is encountered
        self.del_buffer.clear(); // buffer for deletion base series; committed to dd_tag when next match/insertion is encountered
        self.outcomes = PERFECT_MATCH;

        // maintain the same read length as the target strand for best RE-site matching
        // record terminal N bases using the ? operation in the dd_tag
        if aln.tgt_start0 > 0 { 
            seq.push_str('N'.to_string().repeat(aln.tgt_start0).as_str());
            dd_tag.push_str(&format!("?{}", aln.tgt_start0));
            self.outcomes |= HAS_END_CLIP;
        }

        // assemble seq and dd_tag by walking the span of aligned bases
        aln.qry_on_tgt.iter().enumerate().for_each(|(tgt_aln_offset, op_bases)| {
            let tgt_offset = aln.tgt_start0 + tgt_aln_offset;
            let tgt_base  = &strand_pair.this.seq[tgt_offset..=tgt_offset];

            // base identity between strands
            if op_bases == tgt_base {
                self.commit_deletion_run(&mut dd_tag);
                seq.push_str(tgt_base);
                self.ident_len += 1;
                kinetics::push(&strand_pair.this, tgt_offset, false, tx_kinetics);
                qry_offset += 1;
            } else {
                self.commit_identity_run(&mut dd_tag);

                // insertion on qry relative to tgt
                let n_op_bases = op_bases.len();
                if n_op_bases > 1 { 
                    self.commit_deletion_run(&mut dd_tag);
                    // keep the base that the insertion bases were recorded on as is, it is always a match
                    seq.push_str('N'.to_string().repeat(n_op_bases - 1).as_str());
                    dd_tag.push_str(&format!("+{}", op_bases[0..n_op_bases - 1].to_string()));
                    seq.push_str(tgt_base);
                    self.ident_len = 1;
                    qry_offset += n_op_bases;
                    self.outcomes |= HAS_INDEL;
                
                // deletion on qry relative to tgt
                } else if op_bases == "-" { 
                    seq.push('N'); // again, bases were present ("inserted") on one strand but not the other
                    self.del_buffer.push_str(tgt_base);
                    self.outcomes |= HAS_INDEL;
                    // does not advance on qry/prev_strand.seq
                
                // base substitution between strands
                } else { 
                    self.commit_deletion_run(&mut dd_tag);
                    seq.push('N'); // substitution on one strand relative to the other
                    // record the substitution base plus one-base context on either side
                    // note that qry3 is the reverse-complement of tgt3 to match prev_strand ip/pw orientation
                    // thus heteroduplex alignment of ATC (qry = rc(prev_strand.seq)) to ACC (tgt) 
                    // yields tag value *ACCGAT, which parses as:
                    //            ***  <<< the three frames of kinetics data used for assessment
                    //           1 2 3
                    //     5' ---A C C--- tgt    = this_strand.seq
                    //     3' ---T A G--- qry_rc = prev_strand.seq; either strand might be reference
                    //        ---6 5 4---
                    //            ***
                    let (tgt3, tgt_sk) = kinetics::push(
                        &strand_pair.this, tgt_offset, true, tx_kinetics
                    );
                    let (qry3, qry_sk) = kinetics::push(
                        &strand_pair.prev, n_qry_bases - qry_offset - 1, true, tx_kinetics
                    );
                    dd_tag.push_str(&format!("*{}{}", tgt3, qry3));
                    sk_tag.extend(tgt_sk);
                    sk_tag.extend(qry_sk);
                    qry_offset += 1;
                    self.outcomes |= HAS_SUBSTITUTION;
                }
            }
        });

        // finish out the 3' end of the read
        self.commit_identity_run(&mut dd_tag);
        self.commit_deletion_run(&mut dd_tag);
        if aln.tgt_end0 < n_tgt_bases - 1 {
            let n = n_tgt_bases - 1 - aln.tgt_end0;
            seq.push_str('N'.to_string().repeat(n).as_str());
            dd_tag.push_str(format!("?{}", n).as_str());
            self.outcomes |= HAS_END_CLIP;
        }

        // create a two-level QUAL string with Phred 40 and 0 for agreeing and disagreeing bases
        let qual = seq.chars().map(|b| if b == 'N' { 0 } else { 40 }).collect();
        (seq, qual, dd_tag, sk_tag, self.outcomes)
    }

    // commit any buffered identical run to the dd_tag
    fn commit_identity_run(&mut self, dd_tag: &mut String){
        if self.ident_len > 0 {
            dd_tag.push_str(&format!(":{}", self.ident_len));
            self.ident_len = 0;
        }
    }
    // commit any buffered deletion run to the dd_tag
    fn commit_deletion_run(&mut self, dd_tag: &mut String){
        if !self.del_buffer.is_empty() {
            dd_tag.push_str(&format!("+{}", self.del_buffer)); // yes, +, see above
            self.del_buffer.clear();
        } 
    }

}

//! Support for comparing and merging two PacBio strands into a single 
//! duplex consensus read.

// dependencies
use std::collections::HashMap;
use crossbeam::channel::Sender;
use minimap2::{Aligner as Minimap2, Built, Strand};
use faimm::IndexedFasta;
use mdi::pub_key_constants;
use genomex::sequence::{rc_acgt_str, rc_acgtn_str, Aligner, AlignmentStatus};
use crate::tools::basecall_pacbio::channel_kinetics::kinetics;
use super::{StrandPair, KineticsInstance, MAX_READ_LEN};

// constants
pub_key_constants!(
    ALIGNMENT_NOT_FOUND
    ALIGNMENT_TOO_POOR
    ALIGNMENT_TOO_SHORT
);
const DD_CAPACITY: usize = 1000;
const MAX_SHIFT: usize = 15; // allowed Smith-Waterman alignment shift from identity
const MAX_SCORE_DIFF_PER_KB: i32 = 25; // require a minimal match between strands ...
const MAX_SUMMED_CLIPS: usize = 25;    // ... that extends nearly to the read ends
const MIN_REF_ALN_LEN: usize = 100;    // minimum ref alignment length to consider this strand as reference
const PERFECT_MATCH: u8    = 0;

const HAS_END_CLIP: u8     = 1;
const HAS_INDEL: u8        = 2;
const HAS_SUBSTITUTION: u8 = 4;

const MIN_MAPQ: u32 = 40;
const REF_CLIP_OP:         &str = "!"; // marks end clips of  ref_on_this where this.seq did not match reference bases
const DD_CLIP_OP:          &str = "?"; // marks end clips of prev_on_this where this.seq did not match prev.seq bases
const DD_IDENTITY_OP:      &str = ":"; // dd tag operator for homoduplex bases, regardless of reference match
const DD_INSERTION_OP:     &str = "+"; // heteroduplex indels that DO match reference on at least one strand
const DD_DELETION_OP:      &str = "-"; //    INS|DEL identifies the change on the non-reference strand
const DD_INDEL_OP:         &str = "#"; // heteroduplex indels that do NOT match reference on either strand (omitted from seq)
const DD_THIS_MATCH_OP:    &str = ">"; // this.seq base (listed first  in the op value) matches the reference base
const DD_PREV_MATCH_OP:    &str = "<"; // prev.seq base (listed second in the op value) matches the reference base
const DD_MISMATCH_OP:      &str = "*"; // heteroduplex substitutions that do NOT match reference on either strand (committed as N in seq)
const SEQ_MASKED_BASE:     char = 'N'; // base to use when neither strand matches the reference
/// StrandMerger structure for merging two PacBio strands.
pub struct StrandMerger{
    aligner:        Aligner,
    prev_on_this:   Vec<String>,
    ref_on_this:    Vec<String>,
    prev_offset:    usize,
    ident_len:      usize,
    dd_ops:         Vec<(&'static str, String)>,
    is_fusable_op:  HashMap<&'static str, bool>,
    outcomes:       u8,
}
impl StrandMerger {
    /// Create a new StrandMerger structure.
    pub fn new() -> Self {
        StrandMerger { 
            aligner: Aligner::new_fast(
                MAX_READ_LEN + 1,
                MAX_READ_LEN + 1,
                MAX_SHIFT,
            ),
            prev_on_this: Vec::with_capacity(MAX_READ_LEN),
            ref_on_this:  Vec::with_capacity(MAX_READ_LEN),
            prev_offset:  0,
            ident_len:    0,
            dd_ops:       Vec::with_capacity(DD_CAPACITY),
            is_fusable_op: HashMap::from([
                (DD_CLIP_OP,        true),
                (DD_IDENTITY_OP,    false),
                (DD_INSERTION_OP,   true),
                (DD_DELETION_OP,    true),
                (DD_INDEL_OP,       true),
                (DD_THIS_MATCH_OP,  false),
                (DD_PREV_MATCH_OP,  false),
                (DD_MISMATCH_OP,    false),
            ]),
            outcomes:     PERFECT_MATCH,
        }
    }

    /// Align rc(prev.seq) to this.seq to assemble prev_on_this map.
    pub fn set_prev_on_this(
        &mut self,
        strand_pair: &StrandPair,
    ) -> Option<&'static str> { // return Option<failure_reason>

        // align rc(prev.seq) to this.seq 
        let prev = rc_acgt_str(&strand_pair.prev.seq);
        let this = &strand_pair.this.seq;
        let this_len = this.len();
        let mut aln = self.aligner.align(&prev, this, None, true);

        // handle the unxpected outcome of failed alignment between strands
        if aln.status != AlignmentStatus::AlignmentFound {
            return Some(ALIGNMENT_NOT_FOUND);
        }

        // handle poor alignments between strands that might occur from failed strands or other poor ccs
        let max_len = this_len.max(prev.len()) as i32; // the expected perfect alignment score
        let min_allowed_score = max_len - (MAX_SCORE_DIFF_PER_KB * (max_len + 999) / 1000);
        if aln.score < min_allowed_score {
            return Some(ALIGNMENT_TOO_POOR);
        }
        if aln.tgt_end0 - aln.tgt_start0 + 1 < max_len as usize - MAX_SUMMED_CLIPS {
            return Some(ALIGNMENT_TOO_SHORT);
        }

        // fill out the complete prev_on_this map, including terminal clips
        self.prev_on_this.resize(this_len, DD_CLIP_OP.to_string()); // terminal clip operator
        self.prev_on_this[aln.tgt_start0..=aln.tgt_end0].swap_with_slice(aln.qry_on_tgt.as_mut_slice());
        self.prev_offset = aln.qry_start0;
        None
    }

    /// Align this.seq to reference genome to assemble ref_on_this map.
    pub fn set_ref_on_this(
        &mut self,
        strand_pair: &StrandPair,
        minimap2:    &Minimap2<Built>,
        fa:          &IndexedFasta,
    ){

        // perform initial mapping of this.seq to the reference genome
        let this_maps = minimap2.map(
            strand_pair.this.seq.as_bytes(), 
            false, false,  // TODO: perform detailed alignmnent? we will repeat it below
            None, None, None
        ).expect("Minimap2 mapping failure in strand merger");

        // recast the this_to_ref alignment to create ref_on_this, in parallel to structure of prev_on_this
        self.ref_on_this.resize( strand_pair.this.seq.len(), REF_CLIP_OP.to_string());
        for map in this_maps {
            if !map.is_primary || map.mapq < MIN_MAPQ { continue; }
            let this_start = map.query_start  as usize; // half open
            let this_end   = map.query_end    as usize;
            let ref_start  = map.target_start as usize;
            let ref_stop   = map.target_end   as usize;
            if this_end.saturating_sub(this_start) < MIN_REF_ALN_LEN || 
               ref_stop.saturating_sub(ref_start)  < MIN_REF_ALN_LEN { 
                continue; 
            }
            let strand    = map.strand;
            let chrom = map.target_name;
            // let ref_tid    = map.target_id    as usize; // TODO: is this the same as fa tid?
            if let Some(chrom) = chrom {
                if let Some(chrom_index) = fa.fai().tid(&chrom) {
                    let mut ref_seq = fa.view(chrom_index, ref_start, ref_stop)
                        .expect("Failed to extract ref_seq in strand_merger")
                        .to_string()
                        .to_ascii_uppercase();
                    if strand == Strand::Reverse { // this.seq stays the same, we reverse-complement ref_seq to match
                        ref_seq = rc_acgtn_str(&ref_seq);
                    }
                    let mut aln = self.aligner.align(&ref_seq, &strand_pair.this.seq[this_start..this_end], None, true);
                    if aln.status != AlignmentStatus::AlignmentFound {
                        panic!("Failed alignment of faimm 'ref' segment to 'this' strand in strand_merger");
                    }
                    let aln_start = this_start + aln.tgt_start0;
                    let aln_end   = this_start + aln.tgt_end0;
                    self.ref_on_this[aln_start..=aln_end].swap_with_slice(aln.qry_on_tgt.as_mut_slice());
                }
            }
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
        tx_kinetics: &Sender<KineticsInstance>,
    ) -> (String, Vec<u8>, String, Vec<u16>, u8) {

        // initialize merger state
        self.ident_len = 0; 
        self.dd_ops.clear();
        self.outcomes = PERFECT_MATCH;
        let this_len = strand_pair.this.seq.len();
        let prev_len = strand_pair.prev.seq.len();
        let mut seq = String::with_capacity(this_len * 2);
        let mut sk_tag: Vec<u16> = Vec::with_capacity(this_len / 4);

        // loop through this.seq
        // there is one corresponding entry in each of prev_on_this and ref_on_this
        for this_offset in 0..this_len {
            
            // collect the relevant bases and flags
            let this_base = &strand_pair.this.seq[this_offset..=this_offset];
            let prev_bases = &self.prev_on_this[this_offset];
            let this_is_ref = this_base == self.ref_on_this[this_offset];
            let is_heteroduplex = this_base != prev_bases;

            // homoduplex bases that agree between strands
            // commit as is regardless of ref match
            if !is_heteroduplex {
                seq.push_str(this_base);
                self.ident_len += 1;
                if this_is_ref { // only record kinetics benchmarks for bases fully validated by both strands and ref
                    kinetics::push(&strand_pair.this, this_offset, is_heteroduplex, this_is_ref, tx_kinetics);
                }
                self.prev_offset += 1;

            // heteroduplex bases that differ between strands
            // commit the strand that matches the reference, if any
            } else {

                // finish any previous identical run
                if self.ident_len > 0 {
                    self.dd_ops.push((DD_IDENTITY_OP, self.ident_len.to_string()));
                    self.ident_len = 0;
                }
                let prev_is_ref = prev_bases == &self.ref_on_this[this_offset];
                let n_prev_bases = prev_bases.len();

                // short-circuit end clips
                if prev_bases == DD_CLIP_OP {
                    self.dd_ops.push((DD_CLIP_OP, SEQ_MASKED_BASE.to_string()));
                    if this_is_ref {
                        seq.push_str(this_base);
                    } else {
                        seq.push(SEQ_MASKED_BASE);
                    }
                    // prev.offset starts after a left clip
                    self.outcomes |= HAS_END_CLIP;

                // insertion on prev relative to this
                } else if n_prev_bases > 1 {
                    // keep the base that the insertion bases were recorded on as is, it is always a match
                    if this_is_ref {
                        self.dd_ops.push((DD_INSERTION_OP, prev_bases[0..n_prev_bases - 1].to_string()));
                        seq.push_str(this_base);
                    } else if prev_is_ref {
                        self.dd_ops.push((DD_DELETION_OP, prev_bases[0..n_prev_bases - 1].to_string()));
                        seq.push_str(prev_bases);
                    } else {
                        self.dd_ops.push((DD_INDEL_OP, prev_bases[0..n_prev_bases - 1].to_string()));
                        // don't commit, false insertions are more frequent
                        seq.push_str(this_base);
                    }
                    self.ident_len = 1;
                    self.prev_offset += n_prev_bases;
                    self.outcomes |= HAS_INDEL;

                // deletion on prev relative to this
                } else if prev_bases == DD_DELETION_OP { 
                    if this_is_ref {
                        self.dd_ops.push((DD_DELETION_OP, this_base.to_string()));
                        seq.push_str(this_base);
                    } else if prev_is_ref {
                        self.dd_ops.push((DD_INSERTION_OP, this_base.to_string()));
                    } else {
                        self.dd_ops.push((DD_INDEL_OP, this_base.to_string()));
                        // don't commit, false insertions are more frequent
                    }
                    // does not advance on prev_strand.seq
                    self.outcomes |= HAS_INDEL;
                
                // base substitution between strands
                } else { 

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
                    let (this3, this_sk) = kinetics::push(
                        &strand_pair.this, this_offset, is_heteroduplex, this_is_ref, tx_kinetics
                    );
                    let (prev3, prev_sk) = kinetics::push(
                        &strand_pair.prev, prev_len - self.prev_offset - 1, is_heteroduplex, prev_is_ref, tx_kinetics
                    );
                    let dd_op: &str;
                    if this_is_ref {
                        dd_op = DD_THIS_MATCH_OP;
                        seq.push_str(this_base);
                    } else if prev_is_ref {
                        dd_op = DD_PREV_MATCH_OP;
                        seq.push_str(prev_bases);
                    } else {
                        dd_op = DD_MISMATCH_OP;
                        seq.push(SEQ_MASKED_BASE);
                    }
                    self.dd_ops.push((dd_op, format!("{}{}", this3, prev3)));
                    sk_tag.extend(this_sk);
                    sk_tag.extend(prev_sk);
                    self.prev_offset += 1;
                    self.outcomes |= HAS_SUBSTITUTION;
                }
            }
        }

        // finish up any pending identical run
        if self.ident_len > 0 {
            self.dd_ops.push((DD_IDENTITY_OP, self.ident_len.to_string()));
            self.ident_len = 0;
        }

        // collapse runs of the same op to create the final dd tag
        let mut dd_tag = String::with_capacity(DD_CAPACITY);
        let mut dd_iter = self.dd_ops.iter_mut().peekable();
        let mut wrk_op = dd_iter.next().unwrap();
        while let Some(next_op) = dd_iter.peek() {
            if wrk_op.0 == next_op.0 && self.is_fusable_op[wrk_op.0] {
                wrk_op.1.push_str(&next_op.1);
                dd_iter.next();
            } else {
                dd_tag.push_str(wrk_op.0);
                dd_tag.push_str(&wrk_op.1);
                wrk_op = dd_iter.next().unwrap();
            }
        }
        dd_tag.push_str(wrk_op.0);
        dd_tag.push_str(&wrk_op.1);

        // create a two-level QUAL string with:
        //   - Phred 40 for homoduplex or reference-validated bases
        //   - Phred  0 for heteroduplex bases where neither strand matched the reference
        let qual = seq.chars().map(|b| if b == SEQ_MASKED_BASE { 0 } else { 40 }).collect();

        // return the results
        (seq, qual, dd_tag, sk_tag, self.outcomes)
    }

    // // create a map in parallel to qry_on_tgt with reference bases on target
    // // tgt must a slice of bases that matches the span of the cs tag
    // fn set_ref_on_tgt(&mut self, cs: &str, tgt: &str) {
    //     self.ref_on_tgt.clear();
    //     self.op_val.clear();
    //     self.ref_ins.clear();
    //     let mut cs_chars = cs.chars().peekable();
    //     let mut tgt_bases = tgt.chars();
    //     while let Some(op) = cs_chars.next() {
    //         match op {
    //             ':' => { // identical base(s) in ref and tgt
    //                 self.op_val.clear();
    //                 while let Some(&c) = cs_chars.peek() {
    //                     if c.is_ascii_digit() {
    //                         self.op_val.push(cs_chars.next().unwrap());
    //                     } else { break; }
    //                 }
    //                 let base = tgt_bases.next().expect("Target shorter than CS tag");
    //                 if self.ref_ins.is_empty() {
    //                     self.ref_on_tgt.push(base.to_string());
    //                 } else { // prepend pending insertion bases to next base
    //                     self.ref_ins.push(base);
    //                     self.ref_on_tgt.push(self.ref_ins.clone());
    //                     self.ref_ins.clear();
    //                 }
    //                 let n_op_bases: usize = self.op_val.parse().unwrap(); // panic on all parse errors
    //                 if n_op_bases > 1 {
    //                     for _ in 1..n_op_bases {
    //                         let base = tgt_bases.next().expect("Target shorter than CS tag");
    //                         self.ref_on_tgt.push(base.to_string());
    //                     }
    //                 }
    //             },
    //             '*' => { // different base in ref and tgt
    //                 let base = cs_chars.next().unwrap().to_ascii_uppercase();
    //                 if self.ref_ins.is_empty() { // prepend pending insertion bases to next base
    //                     self.ref_on_tgt.push(base.to_string());
    //                 } else {
    //                     self.ref_ins.push(base);
    //                     self.ref_on_tgt.push(self.ref_ins.clone());
    //                     self.ref_ins.clear();
    //                 }
    //                 cs_chars.next(); // consume tgt base
    //                 tgt_bases.next();
    //             },
    //             '+' => { // insertion on tgt relative to ref is a deletion on ref relative to tgt
    //                 while let Some(&c) = cs_chars.peek() {
    //                     if c.is_ascii_alphabetic() {
    //                         self.ref_on_tgt.push("-".to_string());
    //                         cs_chars.next(); // consume bases
    //                         tgt_bases.next().expect("Target shorter than CS tag");
    //                     } else { break; }
    //                 }
    //             },
    //             '-' => { // deletion on tgt relative to ref is a insertion on ref relative to tgt; hold in buffer
    //                 while let Some(&c) = cs_chars.peek() {
    //                     if c.is_ascii_alphabetic() {
    //                         self.ref_ins.push(cs_chars.next().unwrap().to_ascii_uppercase());
    //                     } else { break; }
    //                 }
    //             },
    //             _ => panic!("Unexpected CS tag operation: {}", op),
    //         }
    //     }
    // }

}


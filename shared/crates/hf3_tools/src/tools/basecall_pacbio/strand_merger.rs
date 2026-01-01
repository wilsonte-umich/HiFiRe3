//! Support for comparing and mergin two PacBio strands into a single 
//! duplex consensus read.

// dependencies
use std::collections::HashMap;
use rust_htslib::bam::Record as BamRecord;
use genomex::sequence::Alignment;
use crate::formats::hf_tags::*;
use super::kinetics::Kinetics;
use super::bam;

/// BufferedStrand holds the data for a single PacBio strand while it waits
/// for its matching strand to be encountered for merging.
pub struct BufferedStrand {
    pub usable:   bool,
    pub ec:       f32,
    pub seq:      String,
    pub qual:     Vec<u8>, // TOOD: if buffer size is critical, cache average qual only
    pub ip:       Option<Vec<u8>>, // inter-pulse durations
    pub pw:       Option<Vec<u8>>, // pulse widths
}

/// StrandBuffer structure for caching PacBio strands
/// while waiting for their matching strand.
pub struct StrandBuffer {
    movies: Vec<String>,
    pub strand_buffer: HashMap<usize, BufferedStrand>,  
}
impl StrandBuffer {
    /// Create a new, empty StrandBuffer structure.
    pub fn new() -> Self {
        StrandBuffer {
            movies: Vec::new(),
            strand_buffer: HashMap::new(),
        }
    }
    /// Convert a movie name and ZMW hole number into a unique usize strand buffer key.
    pub fn get_key(&mut self, movie: &str, zmw: usize) -> usize {
        let movie_idx = if let Some(idx) = self.movies.iter().position(|m| m == movie) {
            idx
        } else {
            self.movies.push(movie.to_string());
            self.movies.len() - 1
        };
        movie_idx << 32 | zmw
    }
    /// Push a strand into the buffer for later merging.
    pub fn insert(&mut self, key: usize, usable: bool, ec: f32, strand: &BamRecord) {
        self.strand_buffer.insert(key, BufferedStrand {
            usable,
            ec,
            seq:      strand.seq().as_bytes().iter().map(|&c| c as char).collect(),
            qual:     strand.qual().to_vec(),
            ip:       bam::get_tag_u8_vec_opt(strand, INTER_PULSE_DURATION),
            pw:       bam::get_tag_u8_vec_opt(strand, PULSE_WIDTH),
        });
    }
    /// Attempt to remove a strand from the buffer for merging.
    pub fn remove(&mut self, key: &usize) -> Option<BufferedStrand> {
        self.strand_buffer.remove(&key)
    }
}

/// StrandMerger structure for merging two PacBio strands.
pub struct StrandMerger{
    pub kinetics:   Kinetics,
    pub seq:    String,
    pub dd_tag: String,
    ident_len:  usize,
    del_buffer: String,
    pub qual:   Vec<u8>,
}
impl StrandMerger {

    /// Create a new StrandMerger structure.
    pub fn new() -> Self {
        StrandMerger { 
            kinetics:   Kinetics::new(),
            seq:        String::new(),
            dd_tag:     String::new(),
            ident_len:  0,
            del_buffer: String::new(),
            qual:       Vec::new(),
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
        aln:    &Alignment, // qry in this aligment is the reverse-complemented prev_strand, i.e, rc(qry_rc)
        tgt:    &str,
        qry_rc: &str, // reverse-complement of aln qry/tgt, i.e., the same orientation as the original prev_strand.seq
        tgt_ip: Option<Vec<u8>>,
        tgt_pw: Option<Vec<u8>>,
        qry_ip: Option<Vec<u8>>,
        qry_pw: Option<Vec<u8>>,
    ) {
        // initialize merger state
        self.seq        = String::new(); // the new merged sequence, with Ns at disagreement positions
        self.dd_tag     = String::new(); // the strand differences tag, recording bases at N positions in seq
        self.ident_len  = 0;             // length of a current identical run; committed to dd_tag when a difference is encountered
        self.del_buffer = String::new(); // buffer for deletion base series; committed to dd_tag when next match/insertion is encountered
        let mut qry_offset = aln.qry_start0;
        let n_qry_bases = qry_rc.len();

        // maintain the same read length as the target strand for best RE-site matching
        // record terminal N bases using the ? operation in the dd_tag
        if aln.tgt_start0 > 0 { 
            self.seq.push_str('N'.to_string().repeat(aln.tgt_start0).as_str());
            self.dd_tag.push_str(&format!("?{}", aln.tgt_start0));
        }

        // assemble seq and dd_tag by walking the span of aligned bases
        aln.qry_on_tgt.iter().enumerate().for_each(|(j, qry_bases)| {
            let tgt_offset = aln.tgt_start0 + j;
            let tgt_base  = &tgt[tgt_offset..=tgt_offset];
            if qry_bases == tgt_base { // base identity between strands
                self.commit_deletion_run();
                self.seq.push_str(tgt_base);
                self.ident_len += 1;
                self.kinetics.push(tgt, &tgt_ip, &tgt_pw, tgt_offset, false);
            } else {
                self.commit_identity_run();
                if qry_bases.len() > 1 { // insertion on qry relative to tgt
                    self.commit_deletion_run();
                    // keep the base that the insertion bases were recorded on as is, it is always a match
                    let n_ins_bases = qry_bases.len() - 1; // thus excluding the "carrier" base on tgt
                    self.seq.push_str('N'.to_string().repeat(n_ins_bases).as_str());
                    self.dd_tag.push_str(&format!("+{}", qry_bases[0..n_ins_bases].to_string()));
                    self.seq.push_str(tgt_base);
                    self.ident_len = 1;
                } else if qry_bases == "-" { // deletion on qry relative to tgt
                    self.seq.push('N'); // again, bases were present ("inserted") on one strand but not the other
                    self.del_buffer.push_str(tgt_base);
                } else { // base substitution between strands
                    self.commit_deletion_run();
                    self.seq.push('N'); // substitution on one strand relative to the other
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
                    let tgt3 = self.kinetics.push(
                        tgt, &tgt_ip, &tgt_pw, 
                        tgt_offset, true
                    );
                    let qry3 = self.kinetics.push(
                        qry_rc, &qry_ip, &qry_pw, 
                        n_qry_bases - qry_offset, true
                    );
                    self.dd_tag.push_str(&format!("*{}{}", tgt3, qry3));
                }
            }
            qry_offset += qry_bases.len();
        });

        // finish out the 3' end of the read
        self.commit_identity_run();
        self.commit_deletion_run();
        if aln.tgt_end0 < tgt.len() - 1 {
            let n = tgt.len() - 1 - aln.tgt_end0;
            self.seq.push_str('N'.to_string().repeat(n).as_str());
            self.dd_tag.push_str(format!("?{}", n).as_str());
        }

        // create a two-level QUAL string with Phred 40 and 0 for agreeing and disagreeing bases
        self.qual = self.seq.chars().map(|b| if b == 'N' { 0 } else { 40 }).collect();
    }

    // commit any buffered identical run to the dd_tag
    fn commit_identity_run(&mut self){
        if self.ident_len > 0 {
            self.dd_tag.push_str(&format!(":{}", self.ident_len));
            self.ident_len = 0;
        }
    }
    // commit any buffered deletion run to the dd_tag
    fn commit_deletion_run(&mut self){
        if !self.del_buffer.is_empty() {
            self.dd_tag.push_str(&format!("+{}", self.del_buffer)); // yes, +, see above
            self.del_buffer.clear();
        } 
    }

}

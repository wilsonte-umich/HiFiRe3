/// Handling of PacBio three-strand SNV/indel error correction tags.
/// 
/// Three-strand error correction compares two PacBio strand consensuses 
/// (this, prev) to each other and to the reference genome (ref) to determine 
/// the final basecalling output where:
/// - homoduplex bases with Watson-Crick complementary strands are:
///     - committed as sequenced, regardless of the reference (mis)match
///     - committed with the higher base quality of the two strands
///     - allowed to call SNV and indel variants downstream
/// - heteroduplex bases where one strand matches reference are:
///    - committed as the reference base
///    - committed with the base quality of the strand that matched reference
///    - not relevant to variant calling as they were error-corrected to reference
///    - tracked with kinetics data for analyzing the reason for strand differences
/// - heteroduplex bases where neither strand matches reference, or there is no reference, are:
///    - committed as one or more N bases in SEQ
///        - for unresolved heteroduplex substitutions, a single N is reported
///        - for unresolved heteroduplex indels, the N-track length matches the longer strand
///             e.g., NN is reported if one strand reported two bases where the other reported none
///    - committed with base quality 0 over all reported N bases
///    - not allowed to call SNV and indel variants downstream

// relationship between strand base values, read SEQ and QUAL, and tag operations 
// two (this, prev) or three (this, prev, ref) strands are compared to determine 
//     the basecalling result
// coercion to reference and N/! bases are how duplex error correction is manifest
///    whereas homoduplex bases persist with the ability to call variants
//       homoduplex            heteroduplex
//                        indel      substitution
// this  1  1  1  1     1  1  1  1    1  1  1  1   
// prev  1  1  1  1     2  2  2  2    2  2  2  2
// ref   ?  1  2  ?     ?  1  2  3    ?  1  2  3
// read  1  1  1  1     N  1  2  N    N  1  2  N
// qual  X  X  X  X     !  1  2  !    !  1  2  !  X=strand max, 1|2=reported strand, !=0
// dd:Z  :  =  *  ^     !  +  -  #    ?  >  <  &

// dependencies
use genomex::sequence::Alignment;
use super::*;

// strand_merger outcome flag bits
pub const PERFECT_MATCH: u8                = 0;   // the strand sequences were exactly the same, and matched the reference perfectly
pub const HAS_STRAND_CLIP: u8              = 1;   // bases on prev_strand were clipped when aligning to this_strand
pub const HAS_REF_UNALIGNED: u8            = 2;   // bases on this_strand were not aligned to the reference genome, but others were
pub const HAS_HOMODUPLEX_INDEL: u8         = 4;   // the strands agreed on an indel relative to reference
pub const HAS_STRAND_INDEL: u8             = 8;   // the strands have an indel between them, regardless of its match to reference
pub const HAS_STRAND_INDEL_NEITHER_REF: u8 = 16;  // the strands have an indel between them that does not match reference on either strand
pub const HAS_HOMODUPLEX_SUBS: u8          = 32;  // the strands agreed on a base substitution relative to reference
pub const HAS_STRAND_SUBS: u8              = 64;  // the strands have a base substitution between them, regardless of its match to reference
pub const HAS_STRAND_SUBS_NEITHER_REF: u8  = 128; // the strands have a base substitution between them that does not match reference on either strand

// base to use in merged output SEQ when strands differ and neither strand matches the reference
// homoduplex strands always print bases as sequenced regardless of reference match
// heteroduplex strands print reference bases when one strand matches the reference
pub const SEQ_MASKED_BASE: char = 'N'; 

// DD tag operations
//   prev_on_this clip operations
pub const PREV_CLIP_OP: &str = "~";  // marks end clips of prev_on_this where prev.seq did not match this.seq bases

//   homoduplex operations
pub const HOMODUP_UNKNOWN:  &str = ":"; // homoduplex bases that could not be validated against a good reference alignment (none in the read were)
pub const HOMODUP_REF:      &str = "="; // homoduplex bases that DID     match the reference
pub const HOMODUP_NOT_REF:  &str = "*"; // homoduplex bases that DID NOT match the reference
pub const HOMODUP_NOT_ALN:  &str = "^"; // homoduplex bases that did not align to reference (although others in the read did)

//   heteroduplex indel operations
pub const HETERODUP_INDEL_UNKNOWN:     &str = "!"; // heteroduplex indels that could not be validated against a good reference alignment
pub const HETERODUP_INS_VS_REF:        &str = "+"; // heteroduplex indels that DID match reference on at least one strand
pub const HETERODUP_DEL_VS_REF:        &str = "-"; //    INS|DEL, i.e, +|- identifies the change on the non-reference strand relative to reference
pub const HETERODUP_INDEL_NEITHER_REF: &str = "#"; // heteroduplex indels that DID NOT match reference on either strand

//   heteroduplex base substitution operations
pub const HETERODUP_SUBS_UNKNOWN:      &str = "?"; // heteroduplex substitutions that could not be validated against a good reference alignment
pub const HETERODUP_SUBS_THIS_REF:     &str = ">"; // this.seq base (listed first  in the op value) matched the reference base
pub const HETERODUP_SUBS_PREV_REF:     &str = "<"; // prev.seq base (listed second in the op value) matched the reference base
pub const HETERODUP_SUBS_NEITHER_REF:  &str = "&"; // heteroduplex substitutions that DID NOT match reference on either strand 

// variant calling parameters
const INDEL_FLANK_BASES: usize = 2; // calculate indel base quality including this many bases on either side of the event
const MIN_SNV_INDEL_QUAL: u8 = 27;

/// The types of values found in a `dd:Z:` tag mask.
#[repr(u8)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum DdMaskType {
    Homoduplex = 0,
    EndClipped = 1,
    UnresolvedHeteroduplex = 2,
    CorrectedToReference = 3,
}

// support for SNV consensus building
impl SnvChromWorker {

    /// Convert a `dd:Z:` tag into a `Vec<DdMaskType>` indicating whether each 
    /// read position was either:
    ///     - error corrected to reference 
    ///     - reported as an N base at an unresolved heteroduplex
    /// during basecalling. Homoduplex bases can always call variants. 
    /// EndClipped and UnresolvedHeteroduplex bases will never call variants  
    /// since they were reported as N. CorrectedToReference bases will be 
    /// reported in lower case to guard against the rare situation where the 
    /// reference base was the incorrect choice, i.e., when heteroduplex bases 
    /// were encountered at a true SNP. 
    /// 
    /// For reverse strand alignments, the read mask is reversed to match the read 
    /// SEQ order in the BamRecord.
    pub fn get_dd_mask(
        read: &ReadInstance,
    ) -> Vec<DdMaskType> {
        let read_len = read.qual_bytes.len();
        let mut mask: Vec<DdMaskType> = vec![DdMaskType::Homoduplex; read_len];
        let mut offset0: usize = 0;
        let mut chars = read.dd.chars();
        let mut op = chars.next().unwrap();
        let mut val: String = String::with_capacity(128);
        while let Some(char) = chars.next() {
            if char.is_alphanumeric() {
                val.push(char);
            } else {
                Self::add_to_dd_mask(&mut mask, &mut offset0, op, &val);
                op = char;
                val.clear();
            }
        }
        Self::add_to_dd_mask(&mut mask, &mut offset0, op, &val);
        if read.is_reverse { mask.reverse(); }
        mask
    }

    /// Add one dd tag operation to the read mask.
    fn add_to_dd_mask(
        mask:    &mut Vec<DdMaskType>, 
        offset0: &mut usize,
        op:      char, 
        val:     &str, 
    ) {
        match op {
            // prev_on_this clip operations
            //      two-strand validation of a reference variant is impossible
            '~' => Self::set_dd_mask(mask, offset0, DdMaskType::EndClipped, val.parse::<usize>().unwrap()),
            // homoduplex operations
            //      always allowed to call variants, but only * operations are expected to do so
            //      as alignment outcomes will presumably continue to be the same
            ':' => Self::set_dd_mask(mask, offset0, DdMaskType::Homoduplex, val.parse::<usize>().unwrap()),
            '=' => Self::set_dd_mask(mask, offset0, DdMaskType::Homoduplex, val.parse::<usize>().unwrap()),
            '*' => Self::set_dd_mask(mask, offset0, DdMaskType::Homoduplex, 1), // always come one read base at a time
            '^' => Self::set_dd_mask(mask, offset0, DdMaskType::Homoduplex, val.parse::<usize>().unwrap()),
            // heteroduplex indel operations
            //       ! and # never allowed to call variants since they weren't validated by both read strands
            //       + and - are never expected to call variants as they were error-corrected to reference
            '!' => Self::set_dd_mask(mask, offset0, DdMaskType::UnresolvedHeteroduplex, val.len()), // unknown read and op have same length
            '+' => Self::set_dd_mask(mask, offset0, DdMaskType::CorrectedToReference, 0),      // heterodup insertions relative to ref not included in read
            '-' => Self::set_dd_mask(mask, offset0, DdMaskType::CorrectedToReference,   val.len()), // whereas deletions were committed as ref bases
            '#' => Self::set_dd_mask(mask, offset0, DdMaskType::UnresolvedHeteroduplex, val.len()),
            // heteroduplex base substitution operations
            //       ? and & never allowed to call variants since they weren't validated by both read strands
            //       > and < are never expected to call variants as they were error-corrected to reference
            '?' => Self::set_dd_mask(mask, offset0, DdMaskType::UnresolvedHeteroduplex, 1),
            '>' => Self::set_dd_mask(mask, offset0, DdMaskType::CorrectedToReference,   1),
            '<' => Self::set_dd_mask(mask, offset0, DdMaskType::CorrectedToReference,   1),
            '&' => Self::set_dd_mask(mask, offset0, DdMaskType::UnresolvedHeteroduplex, 1),
            _   => panic!("Unexpected operation in dd tag: {}", op),
        }
    }

    /// Update a block of contiguous read positions in the mask to false as 
    /// needed and increment the position offset.
    fn set_dd_mask(
        mask:    &mut Vec<DdMaskType>,
        offset0: &mut usize, 
        mask_type: DdMaskType, 
        len: usize
    ) {
        if mask_type != DdMaskType::Homoduplex {
            if len > 0 {
                for i in *offset0..(*offset0 + len) {
                    mask[i] = mask_type;
                }
            } else if *offset0 > 0 {
                // heteroduplex insertions relative to ref were corrected to ref
                // thus have no corresponding query base, record mask on previous query base
                mask[*offset0 - 1] = mask_type;
            }
        }
        *offset0 += len;
    }

    /// Process a cs:Z:tag to add a read to a growing fragment variant list. 
    pub fn process_cs_tag(
        &self,
        haplotype:    Haplotype,
        aln:          Option<Alignment>,
        cs_tag:       Option<String>,
        ref_pos0_map: Option<&Vec<ChromPos0>>,
        frag_vars:    &mut FragmentVariants,
        encoding:     &mut AlignmentEncoding,
        read_i:       ReadIndex,
        read:         &ReadInstance,
    ) { 
        let tgt_is_hap = haplotype != Haplotype::Unspecified;
        let (
            mask, 
            qual,
            mut tgt_pos0, // position on either chromosome or haplotype consensus
            mut qry_pos0, // position on the query read
            ref_pos0_map, // mapping of haplotype consensus to chromosome reference
        ) = if tgt_is_hap {(
            Self::get_dd_mask(read),
            read.get_top_strand_qual(),
            aln.as_ref().unwrap().tgt_start0 as u32,
            aln.as_ref().unwrap().qry_start0 as u32,
            ref_pos0_map.unwrap(),
        )} else {(
            Vec::new(),
            Vec::new(),
            read.ref_pos0,
            read.qry_pos0,
            &Vec::new()
        )};

        let mut var_tgt_pos0: Option<SeqPos0> = None;
        let mut n_tgt_bases: u32 = 0;
        let mut alt_bases: UppercaseACGTN = String::with_capacity(128);
        let mut alt_qual: Vec<PhredQual> = Vec::with_capacity(128);
        let mut allowed = true;

        let cs: String;
        let mut chars = if tgt_is_hap {
            cs = cs_tag.unwrap();
            cs.chars()
        } else {
            read.cs.chars()
        };
        let mut op = chars.next().unwrap();
        let mut val: String = String::with_capacity(128);
        
        while let Some(char) = chars.next() {
            if char.is_alphanumeric() {
                val.push(char);
            } else {
                self.handle_cs_op(
                    haplotype, tgt_is_hap, frag_vars, encoding, read_i,
                    &mask, &qual, ref_pos0_map,
                    &mut tgt_pos0, &mut qry_pos0, 
                    &mut var_tgt_pos0, &mut n_tgt_bases, &mut alt_bases, 
                    &mut alt_qual, &mut allowed, 
                    op, &val,
                );
                op = char;
                val.clear();
            }
        }
        self.handle_cs_op(
            haplotype, tgt_is_hap, frag_vars, encoding, read_i,
            &mask, &qual, ref_pos0_map,
            &mut tgt_pos0, &mut qry_pos0, 
            &mut var_tgt_pos0, &mut n_tgt_bases, &mut alt_bases, 
            &mut alt_qual, &mut allowed, 
            op, &val,
        );
    }

    /// Process one cs:Z:tag operation to add to the growing fragment variant 
    /// list. 
    fn handle_cs_op(
        &self,
        haplotype:    Haplotype,
        tgt_is_hap:   bool,
        frag_vars:    &mut FragmentVariants,
        encoding:     &mut AlignmentEncoding,
        read_i:       ReadIndex,
        mask:         &[DdMaskType],
        qual:         &[PhredQual],
        ref_pos0_map: &[ChromPos0],
        tgt_pos0:     &mut SeqPos0,
        qry_pos0:     &mut SeqPos0, 
        var_tgt_pos0: &mut Option<SeqPos0>,
        n_tgt_bases:  &mut u32,
        alt_bases:    &mut UppercaseACGTN,
        alt_qual:     &mut Vec<PhredQual>,
        allowed:      &mut bool,
        op:           char, 
        val:          &str, 
    ) {
        match op {

            // :[0-9]+   Identical sequence length
            ':' => {
                if var_tgt_pos0.is_some() { // commit any preceding variant stretch
                    if *allowed {
                        let variant = Variant::new(
                            var_tgt_pos0.unwrap(),
                            *n_tgt_bases,
                            alt_bases,
                            haplotype
                        );
                        let qual = if tgt_is_hap { 
                            Some({
                                alt_qual.iter().map(|&q| q as f64).sum::<f64>() / 
                                alt_qual.len() as f64
                            } as u8)
                        } else { None };
                        frag_vars.insert(variant, read_i, qual);
                    }
                    *var_tgt_pos0 = None;
                    *n_tgt_bases = 0;
                    alt_bases.clear();
                    alt_qual.clear();
                    *allowed = true;
                }
                let len = val.parse::<u32>().unwrap();
                encoding.add_identity(len);
                *qry_pos0 += len;
                *tgt_pos0 += len;
            },

            // *[acgtn][acgtn]   Substitution: target to query
            '*' => {
                //     S
                // rrrrRrrrr
                // qqqqQqqqq
                //     A
                let alt = val[1..=1].to_ascii_uppercase();
                *allowed &= alt != "N";
                *allowed &= if tgt_is_hap {
                    let ref_pos0 = ref_pos0_map[*tgt_pos0 as usize];
                    !self.simple_repeats.binary_search(ref_pos0, 1)
                } else {
                    !self.simple_repeats.binary_search(*tgt_pos0, 1)
                };
                alt_bases.push_str(&alt);
                let low_qual = if tgt_is_hap {
                    let i0 = *qry_pos0 as usize;
                    alt_qual.push(qual[i0]);
                    *allowed &= mask[i0] != DdMaskType::CorrectedToReference;   
                    qual[i0] <= MIN_SNV_INDEL_QUAL
                } else {
                    false
                };
                encoding.add_substitution(&alt, *allowed, low_qual);
                if var_tgt_pos0.is_none() { *var_tgt_pos0 = Some(*tgt_pos0); }
                *n_tgt_bases += 1;
                *qry_pos0 += 1;
                *tgt_pos0 += 1;
            },

            // +[acgtn]+   Insertion to the target
            '+' => {
                //    *INI     insertions may have heteroduplex bases within homoduplex query run
                // rrrr   Rrrr
                // qqqqQqqqqqq
                //    aA Aa
                let n_ins_bases = val.len();
                let alt = val.to_ascii_uppercase();
                *allowed &= !alt.contains("N");
                *allowed &= if tgt_is_hap {
                    let ref_pos0 = ref_pos0_map[*tgt_pos0 as usize - 1];
                    !self.simple_repeats.binary_search(ref_pos0, 2)
                } else {
                    !self.simple_repeats.binary_search(*tgt_pos0 - 1, 2)
                };
                alt_bases.push_str(&alt);
                let low_qual = if tgt_is_hap {
                    let ins_start0 = *qry_pos0 as usize;
                    let ins_end1 = ins_start0 + n_ins_bases;
                    let qual_left0  = ins_start0.saturating_sub(INDEL_FLANK_BASES);
                    let qual_right1 = (ins_end1 + INDEL_FLANK_BASES).min(qual.len());
                    let q = &qual[qual_left0..qual_right1];
                    alt_qual.extend_from_slice(q);
                    *allowed &= mask[ins_start0] != DdMaskType::CorrectedToReference;
                    let avg_qual = {
                        q.iter().map(|&q| q as f64).sum::<f64>() / 
                        q.len() as f64
                    } as u8;
                    avg_qual < MIN_SNV_INDEL_QUAL
                } else {
                    false
                };
                encoding.add_insertion(*allowed, low_qual);
                if var_tgt_pos0.is_none() { *var_tgt_pos0 = Some(*tgt_pos0 - 1); }
                *qry_pos0 += n_ins_bases as u32;
                // no action on n_tgt_bases or tgt_pos0
            },

            // -[acgtn]+   Deletion from the target
            '-' => {
                //     DDD
                // rrrrRrrrrrr
                // qqqq   Qqqq
                //   aA   Aa
                // heteroduplex indels in read strands are always reported as N bases
                // so do not expect heteroduplex indels to lead to falsely missing bases
                let n_del_bases = val.len() as u32;
                *allowed &= if tgt_is_hap {
                    let ref_pos0 = ref_pos0_map[*tgt_pos0 as usize];
                    !self.simple_repeats.binary_search(ref_pos0, n_del_bases)
                } else {
                    !self.simple_repeats.binary_search(*tgt_pos0, n_del_bases)
                };
                let low_qual = if tgt_is_hap {
                    let qry_after_del0 = *qry_pos0 as usize;
                    let qual_left0  = qry_after_del0.saturating_sub(INDEL_FLANK_BASES);
                    let qual_right1 = (qry_after_del0 + INDEL_FLANK_BASES).min(qual.len());
                    let q = &qual[qual_left0..qual_right1];
                    alt_qual.extend_from_slice(q);
                    *allowed &= mask[qry_after_del0 - 1] != DdMaskType::CorrectedToReference;
                    let avg_qual = {
                        q.iter().map(|&q| q as f64).sum::<f64>() / 
                        q.len() as f64
                    } as u8;
                    avg_qual < MIN_SNV_INDEL_QUAL
                } else {
                    false
                };
                encoding.add_deletion(n_del_bases, *allowed, low_qual);
                if var_tgt_pos0.is_none() { *var_tgt_pos0 = Some(*tgt_pos0); }
                *n_tgt_bases += n_del_bases;
                *tgt_pos0    += n_del_bases;
                // no action on qry_pos0, alt_bases, and N check not applicable
            },
            _   => panic!("Unexpected operation in cs tag: {}", op),
        }
    }

    /// Convert a Smith-Waterman Alignment into the equivalent minimap2 cs tag.
    /// TODO: move this to genomex crate.
    pub fn get_cs_tag(
        // &self,
        aln: &Alignment,
        tgt: &str, // the target sequenced that generated the alignment
    ) -> String { 
        //    M operations carry the query base in the array slot (could be a base mismatch)
        //    I operations carry the inserted base prepended to the NEXT target postion
        //    D operations carry "-" in place of the query base that was deleted relative to target
        let mut cs = String::with_capacity(256);
        let mut del_val: String = String::with_capacity(256);
        let mut identity_len = 0_usize;
        for tgt_i0 in aln.tgt_start0..=aln.tgt_end0 {
            let tgt_base = &tgt[tgt_i0..=tgt_i0];
            let aln_val = aln.qry_on_tgt[tgt_i0 - aln.tgt_start0].as_str();
            if tgt_base == aln_val {
                if del_val.len() > 0 {
                    cs.push_str(&format!("-{}", del_val.to_ascii_lowercase()));
                    del_val.clear();
                }
                identity_len += 1;
            } else {
                if identity_len > 0 { 
                    cs.push_str(&format!(":{identity_len}"));
                    identity_len = 0;
                }
                if aln_val == "-" {
                    del_val.push_str(&tgt_base);
                } else {
                    if del_val.len() > 0 {
                        cs.push_str(&format!("-{}", del_val.to_ascii_lowercase()));
                        del_val.clear();
                    }
                    if aln_val.len() > 1 {
                        let ins_bases = &aln_val[0..aln_val.len() - 1];
                        cs.push_str(&format!("+{}", ins_bases.to_ascii_lowercase()));
                        identity_len = 1;
                    } else {
                        cs.push_str(&format!("*{}{}", tgt_base.to_ascii_lowercase(), aln_val.to_ascii_lowercase()));
                    }
                }
            }

        }
        if identity_len > 0 { cs.push_str(&format!(":{identity_len}")) }
        cs
    }
}

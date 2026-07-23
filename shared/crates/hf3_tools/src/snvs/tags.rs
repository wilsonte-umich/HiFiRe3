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
use crate::snvs::ReferenceVariant;
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

// support for SNV consensus building
impl SnvChromWorker {

    /// Process a cs:Z:tag to add a read to the growing fragment variant list. 
    pub fn fill_fragment_variants(
        &self,
        fvars: &mut FragmentReferenceVariants,
        read_i: ReadIndex,
        aln_start0: u32,
        cs: &str,
    ) { 
        let mut ref_pos0 = aln_start0;
        let mut var_ref_pos0: Option<u32> = None;
        let mut n_ref_bases: u32 = 0;
        let mut alt_bases: String = String::with_capacity(128);
        let mut allowed = true;
        let mut chars = cs.chars();
        let mut op = chars.next().unwrap();
        let mut val: String = String::with_capacity(128);
        while let Some(char) = chars.next() {
            if char.is_alphanumeric() {
                val.push(char);
            } else {
                self.handle_cs_op(
                    fvars, read_i, 
                    &mut ref_pos0, &mut var_ref_pos0, &mut n_ref_bases, 
                    &mut alt_bases, &mut allowed, 
                    op, &val,
                );
                op = char;
                val.clear();
            }
        }
        self.handle_cs_op(
            fvars, read_i, 
            &mut ref_pos0, &mut var_ref_pos0, &mut n_ref_bases, 
            &mut alt_bases, &mut allowed, 
            op, &val,
        );
    }

    /// Process one cs:Z:tag operation to add to the growing fragment variant list. 
    fn handle_cs_op(
        &self,
        fvars: &mut FragmentReferenceVariants,
        read_i: ReadIndex,
        ref_pos0:     &mut u32,
        var_ref_pos0: &mut Option<u32>,
        n_ref_bases:  &mut u32,
        alt_bases:    &mut String,
        allowed:      &mut bool,
        op:  char, 
        val: &str, 
    ) {
        match op {

            // :[0-9]+	Identical sequence length
            ':' => {
                if var_ref_pos0.is_some() { // commit any preceding variant stretch
                    if *allowed {
                        let reference_variant = ReferenceVariant::new(
                            var_ref_pos0.unwrap(),
                            *n_ref_bases,
                            alt_bases,
                        ); 
                        fvars.insert(reference_variant, read_i);
                    }
                    *var_ref_pos0 = None;
                    *n_ref_bases = 0;
                    alt_bases.clear();
                    *allowed = true;
                }
                let len = val.parse::<u32>().unwrap();
                *ref_pos0 += len;
            },

            // *[acgtn][acgtn]	Substitution: ref to query
            '*' => {
                //     S
                // rrrrRrrrr
                // qqqqQqqqq
                //     A
                let alt = val.to_ascii_uppercase();
                *allowed &= alt != "N";
                *allowed &= !self.simple_repeats.binary_search(*ref_pos0, 1);
                alt_bases.push_str(&alt);
                if var_ref_pos0.is_none() { *var_ref_pos0 = Some(*ref_pos0); }
                *n_ref_bases += 1;
                *ref_pos0 += 1;
            },

            // +[acgtn]+	Insertion to the reference
            '+' => {
                //    *INI     insertions may have heteroduplex bases within homoduplex query run
                // rrrr   Rrrr
                // qqqqQqqqqqq
                //    aA Aa
                let alt = val.to_ascii_uppercase();
                *allowed &= !alt.contains("N");
                *allowed &= !self.simple_repeats.binary_search(*ref_pos0 - 1, 2);
                alt_bases.push_str(&alt);
                if var_ref_pos0.is_none() { *var_ref_pos0 = Some(*ref_pos0 - 1); }
                // no action on n_ref_bases or ref_pos0
            },

            // -[acgtn]+	Deletion from the reference
            '-' => {
                //     DDD
                // rrrrRrrrrrr
                // qqqq   Qqqq
                //   aA   Aa
                // heteroduplex indels in read strands are always reported as N bases
                // so do not expect heteroduplex indels to lead to falsely missing bases
                let n_del_bases = val.len() as u32;
                *allowed &= !self.simple_repeats.binary_search(*ref_pos0, n_del_bases);
                if var_ref_pos0.is_none() { *var_ref_pos0 = Some(*ref_pos0); }
                *n_ref_bases += n_del_bases;
                *ref_pos0 += n_del_bases;
                // no action on alt_bases, and N check not applicable
            },
            _   => panic!("Unexpected operation in cs tag: {}", op),
        }
    }


    // /// Process a cs tag into a vector of base values on each reference position. 
    // pub fn set_read_on_ref(
    //     &mut self,
    //     aln_start0:  u32,
    //     site_start0: u32,
    //     cs: &str,
    // ) {
    //     self.read_on_ref.clear(); // always starts at leftmost site position, even if aln doesn't
    //     if aln_start0 > site_start0 {
    //         self.read_on_ref.extend(repeat_n(None, (aln_start0 - site_start0) as usize))
    //     }
    //     let mut ref_pos0 = aln_start0;
    //     let mut chars = cs.chars();
    //     let mut op = chars.next().unwrap();
    //     let mut val: String = String::with_capacity(10);
    //     while let Some(char) = chars.next() {
    //         if char.is_alphanumeric() {
    //             val.push(char);
    //         } else {
    //             self.handle_op(site_start0, &mut ref_pos0, op, &val);
    //             op = char;
    //             val.clear();
    //         }
    //     }
    //     self.handle_op(site_start0, &mut ref_pos0, op, &val);
    // }

    // /// Process one cs:Z:tag operation to add to the growing variant lists.  
    // fn handle_op(
    //     &mut self,
    //     site_start0: u32,
    //     ref_pos0: &mut u32,
    //     op:  char, 
    //     val: &str, 
    // ) {
    //     match op {

    //         // :[0-9]+	Identical sequence length
    //         ':' => {
    //             for _ in 0..val.parse::<usize>().unwrap() {
    //                 if *ref_pos0 >= site_start0 {
    //                     self.read_on_ref.push(Some("=".to_string()));
    //                 }
    //                 *ref_pos0 += 1;
    //             }
    //         },

    //         // *[acgtn][acgtn]	Substitution: ref to query
    //         '*' => {
    //             //     S
    //             // rrrrRrrrr
    //             // qqqqQqqqq
    //             //     A
    //             let alt_base = val.to_ascii_uppercase();
    //             let read_on_ref = if alt_base == "N" ||
    //                 self.simple_repeats.binary_search(*ref_pos0, 1) {
    //                 None
    //             } else {
    //                 Some(alt_base)
    //             };
    //             if *ref_pos0 >= site_start0 {
    //                 self.read_on_ref.push(read_on_ref);
    //             } 
    //             *ref_pos0 += 1;
    //         },

    //         // +[acgtn]+	Insertion to the reference
    //         '+' => {
    //             //    *INI     insertions may have heteroduplex bases within a homoduplex query run
    //             // rrrr   Rrrr
    //             // qqqqQqqqqqq
    //             //    aA Aa
    //             let alt_bases = val.to_ascii_uppercase();
    //             if alt_bases.contains("N") || 
    //                self.simple_repeats.binary_search(*ref_pos0 - 1, 2) {
    //                 // do nothing; ignore disallowed insertions as they do not occupy a reference position
    //             } else if let Some(prev_base) = self.read_on_ref.pop() {
    //                 if let Some(mut prev_base) = prev_base {
    //                     prev_base.push_str(&alt_bases); // record insertion on previous base string
    //                     self.read_on_ref.push(Some(prev_base));
    //                 } else {
    //                     self.read_on_ref.push(Some(alt_bases));
    //                 }
    //             };
    //         },

    //         // -[acgtn]+	Deletion from the reference
    //         '-' => {
    //             //     DDD
    //             // rrrrRrrrrrr
    //             // qqqq   Qqqq
    //             //   aA   Aa
    //             let n_del_bases = val.len();
    //             // heteroduplex indels in read strands are always reported as N bases
    //             // so do not expect heteroduplex indels to lead to falsely missing bases
    //             let read_on_ref = if self.simple_repeats.binary_search(*ref_pos0, n_del_bases as u32) {
    //                 None
    //             } else {
    //                 Some("-".to_string())
    //             };
    //             for _ in 0..n_del_bases {
    //                 if *ref_pos0 >= site_start0 {
    //                     self.read_on_ref.push(read_on_ref.clone());
    //                 }
    //                 *ref_pos0 += 1;
    //             }
    //         },
    //         _   => panic!("Unexpected operation in cs tag: {}", op),
    //     }
    // }

    // /// Process a cs:Z:tag to add a newly encountered read to the growing
    // /// fragment and variant lists. 
    // pub fn process_aln(
    //     &self, 
    //     worker:       &mut SnvChromWorker,
    //     mut ref_pos0: u32,
    //     re_fragment:  &ReFragment,
    //     read_instance: &mut ReadInstance,
    // ) {
    //     let mut chars = self.0.chars();
    //     let mut op = chars.next().unwrap();
    //     let mut val: String = String::with_capacity(10);
    //     let mut qry_pos0: u32 = 0; // alns enforced upstream to have no clips
    //     let mut var_ref_pos0: Option<u32> = None;
    //     let mut n_ref_bases: u32 = 0;
    //     let mut alt_bases: String = String::with_capacity(10);
    //     let mut allowed = true;
    //     while let Some(char) = chars.next() {
    //         if char.is_alphanumeric() {
    //             val.push(char);
    //         } else {
    //             CsTag::process_op(
    //                 worker,
    //                 &mut qry_pos0,  &mut ref_pos0, 
    //                 op, &val, 
    //                 &mut var_ref_pos0, &mut n_ref_bases, &mut alt_bases, &mut allowed,
    //                 re_fragment, read_instance
    //             );
    //             op = char;
    //             val.clear();
    //         }
    //     }
    //     CsTag::process_op(
    //         worker,
    //         &mut qry_pos0,  &mut ref_pos0, 
    //         op, &val, 
    //         &mut var_ref_pos0, &mut n_ref_bases, &mut alt_bases, &mut allowed,
    //         re_fragment, read_instance
    //     );
    // }

    // /// Process one cs:Z:tag operation to add to the growing variant lists.  
    // fn process_op(
    //     worker:       &mut SnvChromWorker,
    //     qry_pos0:     &mut u32, 
    //     ref_pos0:     &mut u32,
    //     op:           char, 
    //     val:          &str, 
    //     var_ref_pos0: &mut Option<u32>,
    //     n_ref_bases:  &mut u32,
    //     alt_bases:    &mut String,
    //     allowed:      &mut bool,
    //     re_fragment:  &ReFragment,
    //     read_instance: &mut ReadInstance,
    // ) {
    //     // :	[0-9]+	Identical sequence length
    //     // *	[acgtn][acgtn]	Substitution: ref to query
    //     // +	[acgtn]+	Insertion to the reference
    //     // -	[acgtn]+	Deletion from the reference
    //     match op {
    //         ':' => {
    //             // commit any prior variant stretch;
    //             if var_ref_pos0.is_some() {
    //                 let reference_variant = ReferenceVariant::new(
    //                     var_ref_pos0.unwrap(), 
    //                     *n_ref_bases, 
    //                     &alt_bases,
    //                 );
    //                 // worker.fragment_variants.insert(
    //                 //     *re_fragment, 
    //                 //     reference_variant.clone()
    //                 // );
    //                 // read_instance.push(reference_variant, *allowed);
    //                 // // worker.variants.increment(
    //                 // //     var_ref_pos0.unwrap(), 
    //                 // //     *n_ref_bases, 
    //                 // //     &alt_bases,
    //                 // //     sample_bit,
    //                 // //     n_passes,
    //                 // //     *allowed,
    //                 // //     qname,
    //                 // // );
    //                 *var_ref_pos0 = None;
    //                 *n_ref_bases = 0;
    //                 alt_bases.clear();
    //                 *allowed = true;
    //             }
    //             let len = val.parse::<u32>().unwrap();
    //             // encoding.add_ref(len);
    //             *qry_pos0 += len;
    //             *ref_pos0 += len;
    //         },
    //         '*' => {
    //             //     S
    //             // rrrrRrrrr
    //             // qqqqQqqqq
    //             //     A
    //             let alt_base = val.chars().nth(1).unwrap().to_ascii_uppercase();
    //             if var_ref_pos0.is_none() { *var_ref_pos0 = Some(*ref_pos0); }
    //             *n_ref_bases += 1;
    //             alt_bases.push(alt_base);
    //             let i0 = *qry_pos0 as usize;
    //             // *allowed &= mask[i0];
    //             *allowed &= !worker.simple_repeats.binary_search(*ref_pos0, 1);
    //             // encoding.add_subs(alt_base, *allowed);
    //             *qry_pos0 += 1;
    //             *ref_pos0 += 1;
    //         },
    //         '+' => {
    //             //    *INI     insertions may have heteroduplex bases within homoduplex query run
    //             // rrrr   Rrrr
    //             // qqqqQqqqqqq
    //             //    aA Aa
    //             let n_ins_bases = val.len();
    //             if var_ref_pos0.is_none() { *var_ref_pos0 = Some(*ref_pos0 - 1); }
    //             alt_bases.push_str(&val.to_ascii_uppercase());
    //             let ins_start0 = *qry_pos0 as usize;
    //             // let ins_end1 = ins_start0 + n_ins_bases;
    //             // *allowed &= mask[ins_start0] && mask[ins_end1 - 1];
    //             *allowed &= !worker.simple_repeats.binary_search(*ref_pos0 - 1, 2);
    //             // encoding.add_ins(*allowed);
    //             *qry_pos0 += n_ins_bases as u32;
    //         },
    //         '-' => {
    //             //     DDD
    //             // rrrrRrrrrrr
    //             // qqqq   Qqqq
    //             //   aA   Aa
    //             let n_del_bases = val.len();
    //             if var_ref_pos0.is_none() { *var_ref_pos0 = Some(*ref_pos0); }
    //             *n_ref_bases += n_del_bases as u32;
    //             let qry_after_del0 = *qry_pos0 as usize;
    //             // *allowed &= mask[qry_after_del0 - 1] && mask[qry_after_del0];
    //             *allowed &= !worker.simple_repeats.binary_search(*ref_pos0, n_del_bases as u32);
    //             // encoding.add_del(n_del_bases as u32, *allowed);
    //             *ref_pos0 += n_del_bases as u32;
    //         },
    //         _   => panic!("Unexpected operation in cs tag: {}", op),
    //     }
    // }
}

// /// DdTag struct helps parse a dd:Z: tag into a mask of read positions that are 
// /// allowed to call SNV and indel variants downstream.
// /// 
// /// For reverse strand alignments, the read mask is reversed to match the read 
// /// SEQ order in the BamRecord.
// pub struct DdTag(pub String);
// impl DdTag {
//     /// Get a `Vec<bool>` indicating whether each read position is allowed to 
//     /// call SNV/indel variants. 
//     /// 
//     /// All read positions can call variants unless masked to false.
//     pub fn get_read_mask(
//         &self, 
//         read_len:   usize, 
//         is_reverse: bool
//     ) -> Vec<bool> {
//         let mut mask: Vec<bool> = vec![true; read_len]; // mask true == allowed
//         let mut offset0: usize = 0;
//         let mut chars = self.0.chars();
//         let mut op = chars.next().unwrap();
//         let mut val: String = String::with_capacity(10);
//         while let Some(char) = chars.next() {
//             if char.is_alphanumeric() {
//                 val.push(char);
//             } else {
//                 DdTag::add_to_mask(&mut mask, &mut offset0, op, &val);
//                 op = char;
//                 val.clear();
//             }
//         }
//         DdTag::add_to_mask(&mut mask, &mut offset0, op, &val);
//         if is_reverse { mask.reverse(); }
//         mask
//     }

//     /// Add one dd tag operation to the read mask.
//     fn add_to_mask(
//         mask:    &mut Vec<bool>, 
//         offset0: &mut usize,
//         op:      char, 
//         val:     &str, 
//     ) {
//         match op {
//             // prev_on_this clip operations
//             //      two-strand validation of a reference variant is impossible
//             '~' => DdTag::set_mask(mask, offset0, false, val.parse::<usize>().unwrap()),
//             // homoduplex operations
//             //      always allowed to call variants, but only * operations are expected to do so
//             //      as alignment outcomes will presumably continue to be the same
//             ':' => DdTag::set_mask(mask, offset0, true, val.parse::<usize>().unwrap()),
//             '=' => DdTag::set_mask(mask, offset0, true, val.parse::<usize>().unwrap()),
//             '*' => DdTag::set_mask(mask, offset0, true, 1), // always come one read base at a time
//             '^' => DdTag::set_mask(mask, offset0, true, val.parse::<usize>().unwrap()),
//             // heteroduplex indel operations
//             //       ! and # never allowed to call variants since they weren't validated by both read strands
//             //       + and - are never expected to call variants as they were error-corrected to reference
//             '!' => DdTag::set_mask(mask, offset0, false, val.len()), // unknown read and op have same length
//             '+' => DdTag::set_mask(mask, offset0, false, 0),    // heterodup insertions relative to ref not included in read
//             '-' => DdTag::set_mask(mask, offset0, false, val.len()), // whereas deletions were committed as ref bases
//             '#' => DdTag::set_mask(mask, offset0, false, val.len()),
//             // heteroduplex base substitution operations
//             //       ? and & never allowed to call variants since they weren't validated by both read strands
//             //       > and < are never expected to call variants as they were error-corrected to reference
//             '?' => DdTag::set_mask(mask, offset0, false, 1),
//             '>' => DdTag::set_mask(mask, offset0, false, 1),
//             '<' => DdTag::set_mask(mask, offset0, false, 1),
//             '&' => DdTag::set_mask(mask, offset0, false, 1),
//             _   => panic!("Unexpected operation in dd tag: {}", op),
//         }
//     }

//     /// Update a block of contiguous read positions in the mask to false as 
//     /// needed and increment the position offset.
//     fn set_mask(
//         mask:    &mut Vec<bool>,
//         offset0: &mut usize, 
//         allowed: bool, 
//         len:     usize
//     ) {
//         if !allowed && len > 0 {
//             for i in *offset0..(*offset0 + len) {
//                 mask[i] = false;
//             }
//         }
//         *offset0 += len;
//     }
// }

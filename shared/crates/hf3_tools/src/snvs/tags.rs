/// Handling of PacBio three-strand SNV/indel error correction tags.
/// 
/// Three-strand error correction compares two PacBio strand
/// consensuses (this, prev) to each other and to the reference genome 
/// (ref) to determine the final basecalling output where:
/// - homoduplex bases with Watson-Crick complementary strands are:
///     - committed as sequenced, regardless of the reference match
///     - allowed to call SNV and indel variants downstream
/// - heteroduplex bases where one strand matches reference are:
///    - committed as reference
///    - not relevant to variant calling as they were error-corrected to reference
///    - tracked with kinetics data for analyzing the reason for strand differences
/// - heteroduplex bases where neither strand matches reference, or there is no good reference, are:
///    - committed as one or more N bases in SEQ
///    - not allowed to call SNV and indel variants downstream

// relationship between strand base values, read SEQ and QUAL, and tag operations 
//       homoduplex      heteroduplex
//                   indel         substitution
// this  1  1  1  1     1  1  1  1    1  1  1  1   two (this, prev) or three strands (this, prev, ref) are compared to determine the basecalling result
// prev  1  1  1  1     2  2  2  2    2  2  2  2
// ref   ?  1  2  ?     ?  1  2  3    ?  1  2  3
// read  1  1  1  1     N  1  2  N    N  1  2  N   
// qual  I  I  I  I     !  I  I  !    !  I  I  !   coercion to reference and N/! bases are the manifestation of duplex base error correction
// dd:Z  :  =  *  ^     !  +  -  #    ?  >  <  &

// dependencies
use crate::snvs::{SnvChromWorker};

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

/// DdTag struct helps parse a dd:Z: tag into a mask of read positions 
/// that are allowed to call SNV and indel variants downstream.
/// 
/// For reverse strand alignments, the read mask is reversed to match 
/// the read SEQ and QUAL order in the BamRecord.
pub struct DdTag(pub String);
impl DdTag {
    /// Get a `Vec<bool>` indicating whether each read position is
    /// allowed to call SNV/indel variants. 
    pub fn get_read_mask(
        &self, 
        read_len:   usize, 
        is_reverse: bool
    ) -> Vec<bool> {
        let mut mask: Vec<bool> = vec![true; read_len]; // all read positions can call variants unless masked to false below
        let mut offset0: usize = 0;
        let mut chars = self.0.chars();
        let mut op = chars.next().unwrap();
        let mut val: String = String::with_capacity(10);
        while let Some(char) = chars.next() {
            if char.is_alphanumeric() {
                val.push(char);
            } else {
                DdTag::add_to_mask(&mut mask, &mut offset0, op, &val);
                op = char;
                val.clear();
            }
        }
        DdTag::add_to_mask(&mut mask, &mut offset0, op, &val);
        if is_reverse { mask.reverse(); }
        mask
    }

    /// Add one dd tag operation to the read mask.
    fn add_to_mask(
        mask:    &mut Vec<bool>, 
        offset0: &mut usize,
        op:      char, 
        val:     &str, 
    ) {
        match op {
            // prev_on_this clip operations
            //      two-strand validation of a reference variant is impossible
            '~' => DdTag::set_mask(mask, offset0, false, val.parse::<usize>().unwrap()),
            // homoduplex operations
            //      always allowed to call variants, but only * operations are expected to do so
            //      as alignment outcomes will presumably continue to be the same
            ':' => DdTag::set_mask(mask, offset0, true, val.parse::<usize>().unwrap()),
            '=' => DdTag::set_mask(mask, offset0, true, val.parse::<usize>().unwrap()),
            '*' => DdTag::set_mask(mask, offset0, true, 1), // always come one read base at a time
            '^' => DdTag::set_mask(mask, offset0, true, val.parse::<usize>().unwrap()),
            // heteroduplex indel operations
            //       ! and # never allowed to call variants since they weren't validated by both read strands
            //       + and - are never expected to call variants as they were error-corrected to reference
            '!' => DdTag::set_mask(mask, offset0, false, val.len()), // unknown read and op have same length
            '+' => DdTag::set_mask(mask, offset0, false, 0),    // heterodup insertions relative to ref not included in read
            '-' => DdTag::set_mask(mask, offset0, false, val.len()), // whereas deletions were committed as ref bases
            '#' => DdTag::set_mask(mask, offset0, false, val.len()),
            // heteroduplex base substitution operations
            //       ? and & never allowed to call variants since they weren't validated by both read strands
            //       > and < are never expected to call variants as they were error-corrected to reference
            '?' => DdTag::set_mask(mask, offset0, false, 1),
            '>' => DdTag::set_mask(mask, offset0, false, 1),
            '<' => DdTag::set_mask(mask, offset0, false, 1),
            '&' => DdTag::set_mask(mask, offset0, false, 1),
            _   => panic!("Unexpected operation in dd tag: {}", op),
        }
    }

    /// Update a block of contiguous read positions in the mask to false
    /// as needed and increment the position offset.
    fn set_mask(
        mask:    &mut Vec<bool>,
        offset0: &mut usize, 
        allowed: bool, 
        len:     usize
    ) {
        if !allowed && len > 0 {
            for i in *offset0..(*offset0 + len) {
                mask[i] = false;
            }
        }
        *offset0 += len;
    }
}

/// CsTag struct helps parse a cs:Z: tag into a read pileup and 
/// allowed variant list.
pub struct CsTag(pub String);
impl CsTag {

    /// Process a cs:Z:tag to add to a pileup and allowed variants list. 
    pub fn process_aln(
        &self, 
        worker:       &mut SnvChromWorker,
        mask:         &[bool],
        mut qry_pos0: u32, 
        mut ref_pos0: u32,
        sample_bit:   u16,
        n_passes:     u8,
    ) {
        let mut chars = self.0.chars();
        let mut op = chars.next().unwrap();
        let mut val: String = String::with_capacity(10);
        let mut var_ref_pos0: Option<u32> = None;
        let mut n_ref_bases: u32 = 0;
        let mut alt_bases: String = String::with_capacity(10);
        let mut allowed = true;
        while let Some(char) = chars.next() {
            if char.is_alphanumeric() {
                val.push(char);
            } else {
                CsTag::process_op(
                    worker, mask,
                    &mut qry_pos0,  &mut ref_pos0, sample_bit, n_passes,
                    op, &val, 
                    &mut var_ref_pos0, &mut n_ref_bases, &mut alt_bases, &mut allowed
                );
                op = char;
                val.clear();
            }
        }
        CsTag::process_op(
            worker, mask,
            &mut qry_pos0,  &mut ref_pos0, sample_bit, n_passes, 
            op, &val, 
            &mut var_ref_pos0, &mut n_ref_bases, &mut alt_bases, &mut allowed
        );
    }

    /// Process one cs:Z:tag operation to add to a pileup and allowed variants list. 
    fn process_op(
        worker:       &mut SnvChromWorker,
        mask:         &[bool],
        qry_pos0:     &mut u32, 
        ref_pos0:     &mut u32,
        sample_bit:   u16,
        n_passes:     u8,
        op:           char, 
        val:          &str, 
        var_ref_pos0: &mut Option<u32>,
        n_ref_bases:  &mut u32,
        alt_bases:    &mut String,
        allowed:      &mut bool,
    ) {
        // :	[0-9]+	Identical sequence length
        // *	[acgtn][acgtn]	Substitution: ref to query
        // +	[acgtn]+	Insertion to the reference
        // -	[acgtn]+	Deletion from the reference
        match op {
            ':' => {
                // commit any prior variant stretch;
                if var_ref_pos0.is_some() {
                    worker.variants.increment(
                        var_ref_pos0.unwrap(), 
                        *n_ref_bases, 
                        &alt_bases,
                        sample_bit,
                        n_passes,
                        *allowed,
                    );
                    *var_ref_pos0 = None;
                    *n_ref_bases = 0;
                    alt_bases.clear();
                    *allowed = true;
                }
                let len = val.parse::<u32>().unwrap();
                worker.pileup.increment_ref(*ref_pos0 as usize, len as usize);
                *qry_pos0 += len;
                *ref_pos0 += len;
            },
            '*' => {
                //     S
                // rrrrRrrrr
                // qqqqQqqqq
                //     A
                let alt_base = val.chars().nth(1).unwrap().to_ascii_uppercase();
                if var_ref_pos0.is_none() { *var_ref_pos0 = Some(*ref_pos0); }
                *n_ref_bases += 1;
                alt_bases.push(alt_base);
                *allowed &= mask[*qry_pos0 as usize];
                worker.pileup.increment_base(*ref_pos0 as usize, alt_base);
                *qry_pos0 += 1;
                *ref_pos0 += 1;
            },
            '+' => {
                //    *INI     insertions may have heteroduplex bases within homoduplex query run
                // rrrr   Rrrr
                // qqqqQqqqqqq
                //     A A
                if var_ref_pos0.is_none() { *var_ref_pos0 = Some(*ref_pos0 - 1); }
                alt_bases.push_str(&val.to_ascii_uppercase());
                *allowed &= mask[*qry_pos0 as usize] && mask[*qry_pos0 as usize + val.len() - 1];
                worker.pileup.increment_ins(*ref_pos0 as usize - 1, *allowed);
                *qry_pos0 += val.len() as u32;
            },
            '-' => {
                //     DDD
                // rrrrRrrrrrr
                // qqqq   Qqqq
                //    A   A
                let len = val.len() as u32;
                if var_ref_pos0.is_none() { *var_ref_pos0 = Some(*ref_pos0); }
                *n_ref_bases += len;
                *allowed &= mask[*qry_pos0 as usize - 1] && mask[*qry_pos0 as usize];
                worker.pileup.increment_del(*ref_pos0 as usize, len as usize, *allowed);
                *ref_pos0 += len;
            },
            _   => panic!("Unexpected operation in cs tag: {}", op),
        }
    }
}

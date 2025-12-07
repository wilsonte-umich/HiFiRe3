//! Reverse-complement functions for DNA sequences and associate quality scores.

/* ------------------------------------------------------------------
DNA bases
------------------------------------------------------------------ */
/// Reverse-complement a sequence containing only ACGT or acgt bases
/// provided as bytes, i.e., `&[u8]`, using fast bit-twiddling.
/// 
/// No check is performed to ensure that the input sequence does not
/// contain invalid bases.
pub fn rc_acgt_u8(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|c| if c & 2 != 0 { c ^ 4 } else { c ^ 21 })
        .collect()
}
/// Reverse-complement a sequence containing only ACGT or acgt bases
/// provided as a string slice, i.e., `&str`. Calls `rc_acgt_u8` internally.
/// 
/// No check is performed to ensure that the input sequence does not
/// contain invalid bases.
pub fn rc_acgt_str(seq: &str) -> String {
    rc_acgt_u8(seq.as_bytes())
        .iter()
        .map(|&c| c as char)
        .collect()
}
/// Reverse-complement a sequence containing only ACGTN or acgtn bases
/// provided as a string slice, i.e., `&str`.
/// 
/// Panic if the input sequence contains invalid bases.
pub fn rc_acgtn_str(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'N' => 'N',
            'a' => 't',
            'c' => 'g',
            'g' => 'c',
            't' => 'a',
            'n' => 'n',
            _   => panic!("rc_acgtn_str encountered an invalid base character in sequence: {}", c),
        })
        .collect()
}

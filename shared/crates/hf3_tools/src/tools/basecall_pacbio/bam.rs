//! Support functions for manipulating unaligned PacBio BAM records.

// dependencies
use rust_htslib::bam::{Record as BamRecord, record::Aux};

/// Get the Aux value of a specified tag from a BAM record.
/// The tags is expected to be present; otherwise this function panics.
pub (super) fn get_tag_value<'a>(strand: &'a BamRecord, tag: &'static str) -> Aux<'a> {
    strand.aux(tag.as_bytes()).unwrap_or_else(|e|{
        panic!("Failed to extract {tag} tag: {}", e)
    })
}
/// Parse the value of a specified 'XX:i:' tag from a BAM record as u8.
/// The tags is expected to be present; otherwise this function panics.
pub (super) fn get_tag_u8(strand: &BamRecord, tag: &'static str) -> u8 {
    match get_tag_value(strand, tag) {
        Aux::U8(v) => v,
        _ => panic!("{tag} tag is not U8"),
    }
}
/// Parse the value of a specified'XX:f:' tag from a BAM record as f32.
/// The tags is expected to be present; otherwise this function panics.
pub (super) fn get_tag_f32(strand: &BamRecord, tag: &'static str) -> f32 {
    match get_tag_value(strand, tag) {
        Aux::Float(v) => v,
        _ => panic!("{tag} tag is not Float"),
    }
}
// /// Parse the value of a specified 'XX:B:C,' tag from a BAM record as Vec<u8>.
// /// The tags is expected to be present; otherwise this function panics.
// pub (super) fn get_tag_u8_vec(strand: &BamRecord, tag: &'static str) -> Vec<u8> {
//     match get_tag_value(strand, tag) {
//         Aux::ArrayU8    (v) => v.iter().collect(),
//         _ => panic!("{tag} tag is not ArrayU8"),
//     }
// }
/// Parse the value of a specified 'XX:B:C,' tag from a BAM record as Option<Vec<u8>>.
/// Return None if the tag is not present.
pub (super) fn get_tag_u8_vec_opt(strand: &BamRecord, tag: &'static str) -> Option<Vec<u8>> {
    if let Some(v) = strand.aux(tag.as_bytes()).ok() {
        return match v {
            Aux::ArrayU8(v) => Some(v.iter().collect()),
            _ => panic!("{tag} tag is not ArrayU8"),
        }
    }
    None
}

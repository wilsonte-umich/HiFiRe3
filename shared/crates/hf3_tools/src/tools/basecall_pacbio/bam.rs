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
/// Get the Aux value of a specified tag from a BAM record.
/// Return None if the tag is not present.
pub (super) fn get_tag_value_opt<'a>(strand: &'a BamRecord, tag: &'static str) -> Option<Aux<'a>> {
    strand.aux(tag.as_bytes()).ok()
}
/// Parse the value of a specified 'XX:i:' tag from a BAM record as u8.
/// The tags is expected to be present; otherwise this function panics.
pub (super) fn get_tag_u8(strand: &BamRecord, tag: &'static str) -> u8 {
    let aux = get_tag_value(strand, tag);
    match aux {
        Aux::I32(v) => v as u8,
        Aux::U8(v)   => v,
        Aux::I16(v) => v as u8,
        Aux::U32(v) => v as u8,
        Aux::I8(v)   => v as u8,
        Aux::U16(v) => v as u8,
        _ => panic!("{tag} tag is not an integer type: {:?}", aux),
    }
}
// /// Parse the value of a specified'XX:f:' tag from a BAM record as f32.
// /// The tags is expected to be present; otherwise this function panics.
// pub (super) fn get_tag_f32(strand: &BamRecord, tag: &'static str) -> f32 {
//     let aux = get_tag_value(strand, tag);
//     match aux {
//         Aux::Float(v) => v,
//         Aux::Double(v) => v as f32,
//         _ => panic!("{tag} tag is not a float type: {:?}", aux),
//     }
// }
/// Parse the value of a specified'XX:f:' tag from a BAM record as f32.
/// Return 0.0 if the tag is not present.
pub (super) fn get_tag_f32_default(strand: &BamRecord, tag: &'static str) -> f32 {
    let aux = get_tag_value_opt(strand, tag);
    match aux {
        Some(Aux::Float(v)) => v,
        Some(Aux::Double(v)) => v as f32,
        Some(_) => panic!("{tag} tag is not a float type: {:?}", aux),
        None => 0.0,
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
    if let Some(aux) = strand.aux(tag.as_bytes()).ok() {
        return match aux {
            Aux::ArrayU8(v) => Some(v.iter().collect()),
            Aux::ArrayI32(v) => Some(v.iter().map(|x| x as u8).collect()),
            _ => panic!("{tag} tag is not ArrayU8: {:?}", aux),
        }
    }
    None
}

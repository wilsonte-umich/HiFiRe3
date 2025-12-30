//! Assess the usability of PacBio CCS strands for duplex "unleaded" basecalling
//! based on the PacBio 'ff:i:' (fail flags) and 'ec:f:' (effective coverage) tags.

// dependencies
use rust_htslib::bam::Record as BamRecord;
use crate::formats::hf_tags::{PACBIO_FAIL, PACBIO_EFF_COVERAGE};
use super::bam;

/// CcsFailBits enumerates the failure reasons for PacBio CCS reads
/// as found in the "ff:i:" tag in ccs bam files (even hifi files)
#[repr(u8)]
enum CcsFailBits {
    // LowAccuracy      = 0x1, // never enforced as is; can rescue low QV reads by strand pairing
    ControlRead      = 0x2,
    // SingleStranded   = 0x4, // not a fail bit here, always set on by-strand ccs reads
    NoCcsRead        = 0x8,
    ChimericAdapter  = 0x10,
    MiscalledAdapter = 0x20,
    CloseAdapter     = 0x40,
}
// 0x1	CCS reads with predicted accuracy below QV 20.
// 0x2	Control CCS reads.
// 0x4	Single-stranded CCS reads.
// 0x8	The median full-length subread from molecules that do not produce a CCS read but have at least one full pass.
// 0x10	CCS reads which are a concatenation of the adapter, with possible short non-adapter sequence in between.
// 0x20	CCS reads with miscalled adapter which is enclosed by a sequence and its reverse complement, either spanning to the end.
// 0x40	CCS reads that have one or more adapters close to either end.

// constants
pub const UNUSABLE_BITS: u8 = 
    CcsFailBits::ControlRead      as u8 |
    CcsFailBits::NoCcsRead        as u8 |
    CcsFailBits::ChimericAdapter  as u8 |
    CcsFailBits::MiscalledAdapter as u8 |
    CcsFailBits::CloseAdapter     as u8;
const MIN_EFFECTIVE_COVERAGE: f32 = 3.0; // TODO: expose as options?
pub (super) const MAX_READ_LEN: usize = 25000;

/// Determine whether a by-strand read is usable for merging based on its PacBio tags.
pub (super) fn is_usable(strand: &BamRecord) -> (bool, f32) {
    let ec = bam::get_tag_f32(strand, PACBIO_EFF_COVERAGE);
    if strand.seq().len() > MAX_READ_LEN { return (false, ec); }
    let ff = bam::get_tag_u8(strand, PACBIO_FAIL);
    ( (ff & UNUSABLE_BITS == 0) && ec >= MIN_EFFECTIVE_COVERAGE, ec )
}

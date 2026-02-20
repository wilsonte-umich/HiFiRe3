//! Assess the usability of PacBio CCS strands for duplex "unleaded" basecalling
//! based on the PacBio 'ff:i:' (fail flags) and 'ec:f:' (effective coverage) tags.

// dependencies
use rust_htslib::bam::Record as BamRecord;
use genomex::bam::tags;
use crate::formats::hf3_tags::{PACBIO_FAIL, PACBIO_EFF_COVERAGE};
use crate::tools::basecall_pacbio::{
    MIN_READ_LEN, 
    MAX_READ_LEN,
    FAIL_BY_OUTCOME, HIFI_BY_OUTCOME, BOTH_BY_OUTCOME,
    FAIL_BY_REASON,  HIFI_BY_REASON,  BOTH_BY_REASON,
};

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
const UNUSABLE_BITS: u8 = 
    CcsFailBits::ControlRead      as u8 |
    CcsFailBits::NoCcsRead        as u8 |
    CcsFailBits::ChimericAdapter  as u8 |
    CcsFailBits::MiscalledAdapter as u8 |
    CcsFailBits::CloseAdapter     as u8;
const MIN_EFFECTIVE_STRAND_COVERAGE: f32 = 3.0; // mergable reads have at least twice this coverage; expose as option?

/// Reasons why a PacBio CCS strand may be considered unusable for merging.
#[repr(u8)]
#[derive(PartialEq, Eq)]
pub enum UnusableReason {
    None           = 0,
    PacBioFlagBits = 1,
    CoverageTooLow = 2,
    ReadTooShort   = 3,
    ReadTooLong    = 4,
}
impl UnusableReason {
    /// Return a string representation of the unusable reason for outcome keying.
    pub fn to_str(&self) -> &'static str {
        match self {
            UnusableReason::None           => "USABLE",
            UnusableReason::PacBioFlagBits => "PACBIO_FAIL_BITS",
            UnusableReason::CoverageTooLow => "COVERAGE_TOO_LOW",
            UnusableReason::ReadTooShort   => "READ_TOO_SHORT",
            UnusableReason::ReadTooLong    => "READ_TOO_LONG",
        }
    }
    pub fn get_reasons(&self, other: &Self) -> &'static str {
        if self == other || other == &UnusableReason::None {
            self.to_str()
        } else if self == &UnusableReason::None {
            other.to_str()
        } else {
            "MIXED_UNUSABLE_REASONS"
        }
    }
}
/// The file source of a strand.
#[repr(u8)]
#[derive(PartialEq, Eq, Copy, Clone)]
pub enum StrandSource {
    Fail = 1, // first since FAIL_BAM_FILE is read first
    Hifi = 2,
    Both = 3, // when one strand came from each file
}
impl StrandSource {
    /// Return the combined source of two strands, which may be the same or different.
    pub fn get_sources(&self, other: &Self) -> Self {
        if self == other {
            *self
        } else {
            StrandSource::Both
        }
    }
    /// Return a string representation of the strand source(s) for outcome keying.
    pub fn to_outcome_key(&self) -> &'static str {
        match self {
            StrandSource::Fail => FAIL_BY_OUTCOME,
            StrandSource::Hifi => HIFI_BY_OUTCOME,
            StrandSource::Both => BOTH_BY_OUTCOME,
        }
    }
    /// Return a string representation of the strand source(s) for reason keying.
    pub fn to_reason_key(&self) -> &'static str {
        match self {
            StrandSource::Fail => FAIL_BY_REASON,
            StrandSource::Hifi => HIFI_BY_REASON,
            StrandSource::Both => BOTH_BY_REASON,
        }
    }
}

/// Determine whether a by-strand read is usable for merging based on its PacBio tags.
pub (super) fn is_usable(strand: &BamRecord) -> (bool, UnusableReason, u8, f32) {
    let ff  = tags::get_tag_u8(strand, PACBIO_FAIL);
    let ec = tags::get_tag_f32_default(strand, PACBIO_EFF_COVERAGE, 0.0);
    let read_len = strand.seq_len();
    let unusable_reason = if ff & UNUSABLE_BITS != 0 {
        UnusableReason::PacBioFlagBits
    } else if ec < MIN_EFFECTIVE_STRAND_COVERAGE {
        UnusableReason::CoverageTooLow
    } else if read_len < MIN_READ_LEN {
        UnusableReason::ReadTooShort
    } else if read_len > MAX_READ_LEN {
        UnusableReason::ReadTooLong // could be masked by CoverageTooLow
    } else {
        UnusableReason::None
    };
    (unusable_reason == UnusableReason::None, unusable_reason, ff, ec)
}

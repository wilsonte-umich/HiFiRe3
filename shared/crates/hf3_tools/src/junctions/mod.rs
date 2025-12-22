//! Module to handle SV junction analysis. 
//! 
//! These implementation are used during both `analyze fragments align`
//! as well as downstream of RE site localization and assignment.

// dependencies

// junction failure flags
#[repr(u8)]
#[derive(PartialEq)]
pub enum JxnFailureFlag {
    None           = 0,
    TraversalDelta = 1, // set by analyze_alignments, independent of RE sites or insert sizing
    Noncanonical   = 2,
    FoldbackInv    = 4,
    OntFollowOn    = 8, 
    HasAdapter     = 16,
    SiteMatch      = 32,  // set by analyze_inserts based on RE site matching and insert sizing
    StemLength     = 64,
    // InsertSize     = 128, // this flag if never set, kept here for legacy reasons
}
const PRE_CHIMERIC_JXN_FLAGS: u8 = 
    JxnFailureFlag::TraversalDelta as u8 |
    JxnFailureFlag::Noncanonical   as u8 |
    JxnFailureFlag::FoldbackInv    as u8;
const CHIMERIC_JXN_FLAGS: u8 = 
    JxnFailureFlag::OntFollowOn as u8 |
    JxnFailureFlag::HasAdapter  as u8 |
    JxnFailureFlag::SiteMatch   as u8 |
    JxnFailureFlag::StemLength  as u8;
impl JxnFailureFlag {
    /// Determine if a junction failure flag indicates a chimeric junction.
    pub fn is_chimeric_jxn_u8(flag: u8) -> bool {
        CHIMERIC_JXN_FLAGS & flag != 0
    }
    // /// Determine if a junction failure flag indicates a chimeric junction.
    // pub fn is_chimeric_jxn(flag: JxnFailureFlag) -> bool {
    //     Self::is_chimeric_jxn_u8(flag as u8)
    // } 
    /// Determine if a junction failed traversal or noncanonical filters,
    /// i.e., upstream of RE site matching and insert sizing.
    pub fn is_pre_chimeric_failure_u8(flag: u8) -> bool {
        PRE_CHIMERIC_JXN_FLAGS & flag != 0
    }
    // /// Determine if a junction failed traversal or noncanonical filters,
    // /// i.e., upstream of RE site matching and insert sizing.
    // pub fn is_pre_chimeric_failure(flag: JxnFailureFlag) -> bool {
    //     Self::is_pre_chimeric_failure_u8(flag as u8)
    // }
}

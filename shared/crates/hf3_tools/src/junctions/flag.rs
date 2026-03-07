// Flags to record the evaluation status of SV junctions relative
// to various chimeric and quality filters.

// A JxnFailureFlag in a bitwise u8 encoding of the reason(s) a
// junction failed filters, if any. Zero indicates the junction
// passed all filters.
#[repr(u8)]
#[derive(PartialEq)]
pub enum JxnFailureFlag {
    None           = 0,
    TraversalDelta = 1, // set by analyze_alignments, independent of RE sites or insert sizing
    Noncanonical   = 2,
    FoldbackInv    = 4,
    LowQualIns     = 8, // either ONT follow-on or PacBio low-quality insertion
    HasAdapter     = 16,
    SiteMatch      = 32,  // set by analyze_inserts based on RE site matching and insert sizing
    StemLength     = 64,
    // InsertSize     = 128, // this flag if never set, kept here for legacy reasons
}

// constants to support matching against multiple junction failure flag bits
const PRE_CHIMERIC_JXN_FLAGS: u8 = 
    JxnFailureFlag::TraversalDelta as u8 |
    JxnFailureFlag::Noncanonical   as u8 |
    JxnFailureFlag::FoldbackInv    as u8;
const CHIMERIC_JXN_FLAGS: u8 = 
    JxnFailureFlag::LowQualIns  as u8 |
    JxnFailureFlag::HasAdapter  as u8 |
    JxnFailureFlag::SiteMatch   as u8 |
    JxnFailureFlag::StemLength  as u8;

// methods for matching against JxnFailureFlag values
impl JxnFailureFlag {

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

    /// Determine if a junction failure flag indicates a chimeric junction.
    pub fn is_chimeric_jxn_u8(flag: u8) -> bool {
        CHIMERIC_JXN_FLAGS & flag != 0
    }

    // /// Determine if a junction failure flag indicates a chimeric junction.
    // pub fn is_chimeric_jxn(flag: JxnFailureFlag) -> bool {
    //     Self::is_chimeric_jxn_u8(flag as u8)
    // } 
}

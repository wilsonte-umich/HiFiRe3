//! Module to handle SV junction analysis. 
//! 
//! These implementation are used during both `analyze fragments align`
//! as well as downstream of RE site localization and assignment.

// dependencies

// junction failure flags
#[repr(usize)]
pub enum JxnFailureFlag {
    None           = 0,
    TraversalDelta = 1,
    Noncanonical   = 2,
    FoldbackInv    = 4,
    OntFollowOn    = 8, 
    HasAdapter     = 16,
    SiteMatch      = 32,
}
pub fn is_chimeric_jxn(flag: JxnFailureFlag) -> bool {
    match flag {
        JxnFailureFlag::OntFollowOn => true,
        JxnFailureFlag::HasAdapter  => true,
        JxnFailureFlag::SiteMatch   => true,
        _ => false,
    }
}


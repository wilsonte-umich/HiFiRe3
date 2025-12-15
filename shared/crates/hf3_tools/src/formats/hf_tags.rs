/// Keep track of HiFiRe3 SAM tags created, or to be retained, at
/// various stages of HiFiRe3 data analysis.
/// 
/// Not all tags are necessarily present on all reads, as some tags
/// are only relevant to specific data types (e.g., ONT vs PacBio).
/// 
/// See file hifire3_sam_tags.csv for extended tag descriptions.

// dependencies
use std::{collections::HashMap, vec};

// constants
// map tag keys to SAM tag prefixes
// Basecalling (these tags are generated entirely outside of HiFiRe3 pipelines)
pub const BASE_MODS: &str               = "ML:B:C:"; // Dorado or PacBio 
pub const BASE_MOD_PROBS: &str          = "MM:Z:";
pub const MEAN_BASE_QUAL: &str          = "qs:f:"; // Dorado
pub const CHANNEL: &str                 = "ch:i:";
pub const POD5_READ_NUMBER: &str        = "rn:i:";
pub const POD5_FILE: &str               = "fn:Z:";
// Trimming
pub const TRIM_LENGTHS: &str            = "tl:Z:"; // hf3_tools trim_ont
// Consensus
pub const PACBIO_CONSENSUS: &str        = "pc:Z:"; // hf3_tools make_pacbio_consensus
// Alignment
pub const FASTP_MERGE: &str             = "fm:Z:"; // fastp in analyze fragments align
pub const ALN_SCORE: &str               = "AS:i:"; // minimap2 in analyze fragments align
pub const MISMATCH_DELETION: &str       = "MD:Z:";
pub const DIFFERENCE_STRING: &str       = "cs:Z:";
pub const DIVERGENCE: &str              = "de:f:";
// AlignmentAnalysis
pub const N_READ_BASES: &str            = "nd:i:"; // hf3_tools analyze_alignments
pub const N_REF_BASES: &str             = "nf:i:";
pub const UNMERGED_NODE5: &str          = "po:Z:";
pub const TARGET_CLASS: &str            = "tc:i:";
pub const IS_ON_TARGET: &str            = "to:i:";
pub const ALN_FAILURE_FLAG: &str        = "af:i:";
pub const BLOCK_N: &str                 = "bn:i:";
pub const JXN_FAILURE_FLAG_TMP: &str    = "jx:i:";
pub const JXN_TYPE: &str                = "jt:i:";
pub const ALN_OFFSET: &str              = "ao:i:";
pub const JXN_BASES: &str               = "jb:Z:";
// InsertAnalysis
pub const SITE_INDEX1_1: &str           = "si:i:"; // hf3_tools analyze_inserts
pub const SITE_POS1_1: &str             = "sp:i:";
pub const SITE_DIST_1: &str             = "sd:i:";
pub const SITE_INDEX1_2: &str           = "xi:i:";
pub const SITE_POS1_2: &str             = "xp:i:";
pub const SITE_DIST_2: &str             = "xd:i:";
pub const SEQ_SITE_INDEX1_2: &str       = "qi:i:";
pub const SEQ_SITE_POS1_2: &str         = "qp:i:";
pub const IS_END_TO_END_READ: &str      = "re:i:";
pub const IS_END_TO_END_INSERT: &str    = "ie:i:";
pub const NODE_5: &str                  = "df:i:";
pub const NODE_3: &str                  = "dl:i:";
pub const INSERT_SIZE: &str             = "iz:i:";
pub const IS_ALLOWED_SIZE: &str         = "az:i:";
pub const STEM5_LENGTH: &str            = "ul:i:";
pub const STEM3_LENGTH: &str            = "vl:i:";
pub const PASSED_STEM5: &str            = "up:i:";
pub const PASSED_STEM3: &str            = "vp:i:";
pub const JXN_FAILURE_FLAG: &str        = "jf:i:";
pub const READ_HAS_JXN: &str            = "rj:i:";

/// StageTags enumerate the analysis states at which SAM tags are
/// added, or may be retained, in HiFiRe3 data processing stages.
pub enum StageTags {
    BaseCalling,
    Trimming,
    Consensus,
    Alignment,
    AlignmentAnalysis,
    InsertAnalysis,
}
impl StageTags {
    /// Return a vector of SAM tags that may added by a specified analysis
    /// stage, depending on the type of read data being analyzed.
    pub fn tag_added_by_stage(self) -> Vec<&'static str> {
        match self {
            Self::BaseCalling => vec![
                BASE_MODS,
                BASE_MOD_PROBS,
                MEAN_BASE_QUAL,
                CHANNEL,
                POD5_READ_NUMBER,
                POD5_FILE,
            ],
            Self::Trimming => vec![
                TRIM_LENGTHS,
            ],
            Self::Consensus => vec![
                PACBIO_CONSENSUS,
            ],
            Self::Alignment => vec![
                FASTP_MERGE,
                ALN_SCORE,
                MISMATCH_DELETION,
                DIFFERENCE_STRING,
                DIVERGENCE,
            ],
            Self::AlignmentAnalysis => vec![
                N_READ_BASES,
                N_REF_BASES,
                UNMERGED_NODE5,
                TARGET_CLASS,
                IS_ON_TARGET,
                ALN_FAILURE_FLAG,
                BLOCK_N,
                JXN_FAILURE_FLAG_TMP,
                JXN_TYPE,
                ALN_OFFSET,
                JXN_BASES,
            ],
            Self::InsertAnalysis => vec![
                SITE_INDEX1_1,
                SITE_POS1_1,
                SITE_DIST_1,
                SITE_INDEX1_2,
                SITE_POS1_2,
                SITE_DIST_2,
                SEQ_SITE_INDEX1_2,
                SEQ_SITE_POS1_2,
                IS_END_TO_END_READ,
                IS_END_TO_END_INSERT,
                NODE_5,
                NODE_3,
                INSERT_SIZE,
                IS_ALLOWED_SIZE,
                STEM5_LENGTH,
                STEM3_LENGTH,
                PASSED_STEM5,
                PASSED_STEM3,
                JXN_FAILURE_FLAG,
                READ_HAS_JXN,
            ],
        }
    }

    /// Return a vector of SAM tags that may present after a specified analysis 
    /// stage has completed, depending on the type of read data being analyzed.
    pub fn tags_after_stage(self) -> Vec<&'static str> {
        match self {
            Self::BaseCalling => Self::BaseCalling.tag_added_by_stage(),
            Self::Trimming => {
                let mut tags = Self::BaseCalling.tag_added_by_stage();
                tags.extend(Self::Trimming.tag_added_by_stage());
                tags
            },
            Self::Consensus => {
                let mut tags = Self::BaseCalling.tag_added_by_stage();
                tags.extend(Self::Consensus.tag_added_by_stage());
                tags
            },
            Self::Alignment => {
                let mut tags = Self::Trimming.tags_after_stage();
                tags.extend(Self::Consensus.tag_added_by_stage()); // Trimming and Consensus tags won't both be present
                tags.extend(Self::Alignment.tag_added_by_stage());
                tags
            },
            Self::AlignmentAnalysis => {
                let mut tags = Self::Alignment.tags_after_stage();
                tags.extend(Self::AlignmentAnalysis.tag_added_by_stage());
                tags
            },
            Self::InsertAnalysis => {
                let mut tags = Self::AlignmentAnalysis.tags_after_stage();
                tags.extend(Self::InsertAnalysis.tag_added_by_stage());
                tags
            },
        }
    }
}

// pub struct TagStore {
//     pub tag_map: HashMap<&'static str, bool>,
// }
// impl TagStore {

//     /// Initialize a TagStore for storing and updating tag values prior to final printing.
//     pub fn new(state: TagState) -> TagStore {
//         let tags = state.tag_at_state();
//         let mut tag_map: HashMap<&'static str, bool> = HashMap::new();
//         for tag in tags {
//             tag_map.insert(tag, true);
//         }
//         TagStore{
//             tag_map,
//         }
//     }

//     /// Check if a tag is to be retained.
//     pub fn is_retained(&self, tag: &str) -> bool {
//         self.tag_map.contains_key(tag)
//     }
// }

/// Keep track of HiFiRe3 SAM tags created, or to be retained, at
/// various stages of HiFiRe3 data analysis.
/// 
/// Not all tags are necessarily present on all reads, as some tags
/// are only relevant to specific data types (e.g., ONT vs. PacBio)
/// or read configurations.

/* -----------------------------------------------------------------
map tag keys to SAM tag/aux prefixes
----------------------------------------------------------------- */
// Basecalling (these tags are generated entirely outside of/before HiFiRe3 pipelines)
pub const BASE_MODS: &str               = "ML:B:C,"; // Dorado or PacBio SAM modification tags; only if available
pub const BASE_MOD_PROBS: &str          = "MM:Z:";
pub const CHANNEL: &str                 = "ch:i:";   // Dorado tags; only on ONT reads
pub const POD5_READ_NUMBER: &str        = "rn:i:";   // POD5 lookup tags; only in unaligned BAM, dropped in aligned BAM
pub const POD5_FILE: &str               = "fn:Z:";
pub const PACBIO_FAIL: &str             = "ff:i:";   // PacBio tags; only on PacBio reads
pub const INTER_PULSE_DURATION: &str    = "ip:B:C,"; // kinetics tags as found on by-strand CCS reads
pub const PULSE_WIDTH: &str             = "pw:B:C,";
pub const PACBIO_EFF_COVERAGE: &str     = "ec:f:";   // effective coverage, modified by summing during consensus building

// Trimming, added by hf3_tools trim_ont or reformat_ont
pub const TRIM_LENGTHS: &str            = "tl:Z:";   // serialized encoding of 5',3' ONT trim lengths; only on ONT reads

// Consensus, added by hf3_tools basecall_pacbio
pub const STRAND_DIFFERENCE_TYPES: &str = "dt:i:";   // bit-encoded flag of the types of duplex strand differences; absent if duplex basecalling not performed
pub const STRAND_DIFFERENCES: &str      = "dd:Z:";   // cs-like encoding of strand differences, includes reference-resolved heteroduplex (+, -, >, <)
pub const SUBSTITUTION_KINETICS: &str   = "sk:B:S,"; // kinetics values surrounding base substitution differences

// Alignment, added by fastp, minimap2, or the Ultima Aligner
pub const FASTP_MERGE: &str             = "fm:Z:";   // serialized encoding of read1,read2 retained bases; only on merged read pairs
pub const ALN_SCORE: &str               = "AS:i:";   // minimap2 tags 
pub const DIVERGENCE: &str              = "de:f:";
pub const DIFFERENCE_STRING: &str       = "cs:Z:";   // only retained when needed for variant calling; dropped on perfect reference matches

// AlignmentAnalysis, added by hf3_tools analyze_alignments
pub const PAIRED_OUTER_NODE: &str       = "po:Z:";   // signed 64-bit encoding of the paired read 5' outer node; only on unmerged read pairs
pub const TARGET_MATCH: &str            = "tm:i:";   // bit-encoded target region metadata as region_i1 << 3 | paired_state; only on targeted libraries
pub const READ_IS_OFF_TARGET: &str      = "to:A:";   // flag indicating that the read is considered to be off-target; on-target/untargeted if absent
pub const ALN_FAILURE_FLAG: &str        = "af:i:";   // bit-encoded flag based on MAQ and divergence; zero/passed if absent
pub const BLOCK_N: &str                 = "bn:i:";   // number/index of the traversal block of this alignment; block 1 if absent
pub const JUNCTION: &str                = "jx:Z:";   // serialized encoding of junction metadata; only on alignments followed by another alignment
pub const JXN_FAILURE_FLAG_INIT: &str   = "ji:i:";   // bit-encoded flag based on all junction filters to this point; zero/passed if absent

// InsertAnalysis, added by hf3_tools analyze_inserts
pub const READ_FAILURE_FLAG: &str       = "rf:i:";   // bit-encoded flag why a read was rejected (unmapped, off-target, excessive outer clip, no outer site); zero/passed if absent
pub const JXN_FAILURE_FLAG: &str        = "jf:i:";   // bit-encoded flag based on all junction filters; zero/passed if absent
pub const IS_END_TO_END: &str           = "ee:i:";   // bit-encoded flag 000000IR whether the parent read (R) and insert (I) is end-to-end; on all on-target reads
pub const OUTER_NODES: &str             = "on:Z:";   // paired signed 64-bit encoding of 5',3' outer nodes; on all on-target reads
pub const CLOSEST_SITES: &str           = "sc:B:i,"; // distances to closest 5',3',projected sites, each as signed-index1,site_pos1,signed-dist; only if matching sites
pub const INSERT_SIZE: &str             = "iz:i:";   // signed integer insert size; positive = passed 1N to 2N filter; only on sized libraries
pub const STEM_LENGTH5: &str            = "ul:i:";   // signed integer 5' stem length; positive = passed <1N filter; only on sized libraries
pub const STEM_LENGTH3: &str            = "vl:i:";   // signed integer 3' stem length; positive = passed <1N filter; only on sized libraries
pub const READ_HAS_PASSED_JXN: &str     = "rj:A:";   // flag that the parent read has at least one passed junction; false if absent

// CompareS(N)Vs, added by split_by_chrom_s(n)vs to support sample tracking
pub const SAMPLE_BIT: &str              = "sb:i:";   // bit-encoded sample index, to track the source sample vs. a known sample list

/* -----------------------------------------------------------------
enumerate the tags added or retained at each analysis stage
----------------------------------------------------------------- */
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
                CHANNEL,
                POD5_READ_NUMBER, // POD5 lookup tags only in unaligned BAM
                POD5_FILE,
                PACBIO_FAIL,
                INTER_PULSE_DURATION,
                PULSE_WIDTH,
                PACBIO_EFF_COVERAGE,
            ],
            Self::Trimming => vec![
                TRIM_LENGTHS,
            ],
            Self::Consensus => vec![
                // PACBIO_EFF_COVERAGE,
                STRAND_DIFFERENCE_TYPES,
                STRAND_DIFFERENCES,
                SUBSTITUTION_KINETICS,
            ],
            Self::Alignment => vec![
                FASTP_MERGE,
                ALN_SCORE,
                DIVERGENCE,
                DIFFERENCE_STRING,
            ],
            Self::AlignmentAnalysis => vec![
                PAIRED_OUTER_NODE,
                TARGET_MATCH,
                READ_IS_OFF_TARGET,
                ALN_FAILURE_FLAG,
                BLOCK_N,
                JUNCTION,
                JXN_FAILURE_FLAG_INIT,
            ],
            Self::InsertAnalysis => vec![
                READ_FAILURE_FLAG,
                JXN_FAILURE_FLAG,
                IS_END_TO_END,
                OUTER_NODES,
                CLOSEST_SITES,
                INSERT_SIZE,
                STEM_LENGTH5,
                STEM_LENGTH3,
                READ_HAS_PASSED_JXN,
            ],
        }
    }

    /// Return a vector of SAM tags that may present and are to be retained 
    /// after a specified analysis stage has completed, depending on the 
    /// type of read data being analyzed.
    pub fn tags_after_stage(self) -> Vec<&'static str> {
        match self {
            Self::BaseCalling => vec![
                BASE_MODS,
                BASE_MOD_PROBS,
                CHANNEL,             // drop POD5 lookup tags from alignment BAM
                PACBIO_FAIL,         // ip and pw tags dropped by basecall_pacbio in favor of dd, sk and dt tags
                PACBIO_EFF_COVERAGE, // ff and ec tags updated by basecall_pacbio to reflect both merged strands
            ],
            // retain all pre-alignment tags (BaseCalling, Trimming, Consensus) until Alignment stage
            Self::Trimming => { 
                let mut tags = Self::BaseCalling.tag_added_by_stage(); // NOT tags_after_stage
                tags.extend(Self::Trimming.tag_added_by_stage());
                tags
            },
            Self::Consensus => {
                let mut tags = Self::BaseCalling.tag_added_by_stage(); // NOT tags_after_stage
                tags.extend(Self::Consensus.tag_added_by_stage());
                tags
            },
            // drop some pre-alignment tags at Alignment stage
            Self::Alignment => {
                let mut tags = Self::BaseCalling.tags_after_stage();
                tags.extend(Self::Trimming.tag_added_by_stage());
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

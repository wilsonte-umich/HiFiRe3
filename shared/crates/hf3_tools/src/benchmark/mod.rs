//! Functions to support benchmarking against VCF files. 

// modules
mod re_site;
mod vcf;
mod sv_matching;

// re-exports
pub use re_site::*;
pub use vcf::*;
pub use sv_matching::*;

// constants
pub const MAX_CONSENSUS_SEQS: usize = 25; // to maintain efficiency of consensus building
pub const SVINSSEQ: &str = "SVINSSEQ";
pub const SVINSSEQ_TAG: &[u8] = b"SVINSSEQ";
pub const SCAFFOLD_FLANK_LEN: u32 = 15_000; // pad bases around junctions to size 6-base RE fragments
pub const INSIDE_REGION: &str = "I";

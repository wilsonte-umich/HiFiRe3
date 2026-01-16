
// modules

// ONT read processing
pub mod trim_ont;
pub mod reformat_ont;

// PacBio by-strand ead processing
pub mod basecall_pacbio;

// read alignment analysis
pub mod analyze_alignments;
pub mod analyze_inserts;

// structural variant analysis
pub mod split_by_chrom;
pub mod index_fragments_by_chrom;


// modules

// ONT read processing
pub mod trim_ont;     // the main trimming tool for HiFiRe3 ONT basecalling
pub mod reformat_ont; // alternative to reformat legacy basecall and trim files

// PacBio by-strand read merging
pub mod basecall_pacbio;

// read alignment analysis
pub mod analyze_alignments; // applied to the minimap2 alignment stream
pub mod analyze_inserts;    // applied to all alignments after site localization

// structural variant analysis
pub mod split_bam_by_chrom; // create temporary chrom-level BAM files
pub mod analyze_svs;        // analyze chrom-level BAM files then aggregate SV data

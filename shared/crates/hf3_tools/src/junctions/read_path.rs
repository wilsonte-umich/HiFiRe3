//! Structure and methods for propagating SV read path information
//! to files used by the R Shiny app.

// dependencies
use std::error::Error;
use std::cmp::Ordering;
use std::str::from_utf8_unchecked;
use rayon::prelude::*;
use rust_htslib::{bam::{HeaderView, record::Record as BamRecord}};
use mdi::{InputCsv, OutputCsv};
use genomex::bam::{tags, cigar};
use crate::formats::hf3_tags::*;

// constants
const COMMA: &str = ",";
const PHRED_OFFSET: u8 = 33;

/// A SvReadPath holds detailed information on the path of alignments
/// through a read with a committed junction, i.e., with multiple alignment 
/// segments that passed traversal delta.
#[derive(serde::Serialize, serde::Deserialize)]
pub struct SvReadPath {
    /* ------------------------------------------- */
    // read identifier and read-level metadata
    #[serde(rename(serialize = "#qname", deserialize = "qname"))]
    qname:          String,
    read_len:       u32, // length of SEQ
    insert_size:    i32, // INSERT_SIZE tag value, or read_len if insert_size not available
    has_passed_jxn: u8,  // READ_HAS_PASSED_JXN tag value
    /* ------------------------------------------- */
    // concatenated alignment-level metadata
    // use String up front since SvReadPath only exists to be written to file
    chroms:            String, // chromosome names of the alignment segments
    pos1s:             String, // 1-based leftmost positions of the alignment segments
    strand0s:          String, // strands of the alignment segments as 0=top|1=bottom
    n_ref_bases:       String, // number of reference bases spanned
    qry_start0s:       String, // query start positions
    qry_end1s:         String, // query end   positions
    block_ns:          String, // BLOCK_N tag values
    mapqs:             String, // MAPQ values
    divergences:       String, // DIVERGENCE tag values
    aln_failure_flags: String, // ALN_FAILURE_FLAG tag values
    jxn_failure_flags: String, // JXN_FAILURE_FLAG tag values (0 or last alignment where not applicable)
    cigars:            String, // CIGAR strings
    /* ------------------------------------------- */
    // full read SEQ and QUAL as reported by the first alignment of the read
    // reverse complement of the read if seq_strand0 == 1
    seq_strand0: u8, // strand of seq as 0=top|1=bottom; compare strand0s[i] to seq_strand0 to match seq to cigars[i]
    seq:         String, 
    qual:        String,
}
impl SvReadPath {

    /// Create a new pre-allocated SvReadPath.
    pub fn from_alns(
        alns:    &[BamRecord],
        header:  &HeaderView,
    ) -> Self{
        let n_alns = alns.len();
        let read_len = alns[0].seq_len() as u32;

        // initialize the SvReadPath with read-level values
        let mut read_path = unsafe { SvReadPath{
            qname:            from_utf8_unchecked(alns[0].qname()).to_string(),
            read_len:         read_len,
            insert_size:      tags::get_tag_i32_default( &alns[0], INSERT_SIZE, read_len as i32),
            has_passed_jxn:   tags::get_tag_bool_default(&alns[0], READ_HAS_PASSED_JXN, false) as u8,
            /* ------------------------------------------- */
            chroms:           String::with_capacity(n_alns * 11), // try to avoid reallocations
            pos1s:            String::with_capacity(n_alns * 10),
            strand0s:         String::with_capacity(n_alns * 2),
            n_ref_bases:      String::with_capacity(n_alns * 7),
            qry_start0s:      String::with_capacity(n_alns * 10),
            qry_end1s:        String::with_capacity(n_alns * 10),
            block_ns:         String::with_capacity(n_alns * 3),
            mapqs:            String::with_capacity(n_alns * 3),
            divergences:      String::with_capacity(n_alns * 20),
            aln_failure_flags:String::with_capacity(n_alns * 4),
            jxn_failure_flags:String::with_capacity(n_alns * 4),
            cigars:           String::with_capacity(n_alns * 500),
            /* ------------------------------------------- */
            seq_strand0:      alns[0].is_reverse() as u8,
            seq:              from_utf8_unchecked(&alns[0].seq().as_bytes()).to_string(),
            qual:             alns[0].qual().iter().map(|q| (q + PHRED_OFFSET) as char).collect(),
        } };

        // iterate to fill in concatenated alignment-level values
        // all values, including single and last values, end with a comma
        unsafe { alns.iter().for_each(|aln| {
            read_path.chroms           .push_str(&(from_utf8_unchecked(header.tid2name(aln.tid() as u32)).to_string() + COMMA));
            read_path.pos1s            .push_str(&((aln.pos() + 1).to_string() + COMMA));
            read_path.strand0s         .push_str(&((aln.is_reverse() as u8).to_string() + COMMA));
            read_path.n_ref_bases      .push_str(&((aln.cigar().end_pos() - aln.pos()).to_string() + COMMA));
            read_path.qry_start0s      .push_str(&(cigar::get_query_start0(aln).to_string() + COMMA));
            read_path.qry_end1s        .push_str(&(cigar::get_query_end1(aln, read_len).to_string() + COMMA));
            read_path.block_ns         .push_str(&(tags::get_tag_u8_default(aln, BLOCK_N, 1).to_string() + COMMA));
            read_path.mapqs            .push_str(&(aln.mapq().to_string() + COMMA));
            read_path.divergences      .push_str(&(tags::get_tag_f32(aln, DIVERGENCE).to_string() + COMMA));
            read_path.aln_failure_flags.push_str(&(tags::get_tag_u8_default(aln, ALN_FAILURE_FLAG, 0).to_string() + COMMA));
            read_path.jxn_failure_flags.push_str(&(tags::get_tag_u8_default(aln, JXN_FAILURE_FLAG, 0).to_string() + COMMA));
            read_path.cigars           .push_str(&(aln.cigar().to_string() + COMMA));
        }); }

        // return the filled SvReadPath
        read_path
    }

    /// Write a vector of AlignmentSegments to a file, one row per
    /// unique sorted segment with a count.
    pub fn write_sorted(
        mut read_paths: Vec<SvReadPath>, 
        filepath: &str,
        ncpu:     u32,
    ) -> Result<(), Box<dyn Error>> {
        read_paths.par_sort_unstable();
        let writer = OutputCsv::open(filepath, Some(ncpu));
        writer.serialize_all(&read_paths);
        Ok(())
    }  

    /// Merge two or more read path files.
    pub fn merge_and_write_sorted(
        merge_input_dirs: &[&str],
        filepath: &str,
        ncpu:     u32,
    ) -> Result<(), Box<dyn Error>> {
        let mut read_paths: Vec<SvReadPath> = Vec::new();
        for dir in merge_input_dirs {
            let mut reader = InputCsv::open_file_from_glob(
                dir, "read_paths.txt.bgz", 
                b'\t', true
            )?;
            for result in reader.deserialize(){
                let read_path = result?;
                read_paths.push(read_path);
            }
        }
        SvReadPath::write_sorted(read_paths, filepath, ncpu)
    }
}

// support sorting of read paths by qname and read length
impl PartialEq for SvReadPath {
    fn eq(&self, other: &Self) -> bool {
        self.qname == other.qname
    }
}
impl Eq for SvReadPath {}
impl PartialOrd for SvReadPath {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for SvReadPath {
    fn cmp(&self, other: &Self) -> Ordering {
        self.qname.cmp(&other.qname)
    }
}

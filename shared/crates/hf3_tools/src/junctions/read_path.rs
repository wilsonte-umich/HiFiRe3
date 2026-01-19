//! Structure and methods for propagating SV read path information
//! to files used by the R Shiny app.

// dependencies
use std::error::Error;
use std::io::Write;
use rust_htslib::bgzf::Writer as BgzWriter;
use std::cmp::Ordering;
use std::str::from_utf8_unchecked;
use rust_htslib::bam::{record::Record as BamRecord, HeaderView};
use rayon::prelude::*;
use genomex::bam::{tags, cigar, junction::pack_signed_node_aln};
use genomex::genome::Chroms;
use crate::formats::hf3_tags::*;

// constants
const COMMA: &str = ",";

/// A SvReadPath holds detailed information on the path of alignments
/// through a read with a junction, i.e., with multiple alignment segments.
pub struct SvReadPath {
    /* ------------------------------------------- */
    // read identifier and read-level metadata
    qname:          String,
    read_len:       u32,  // length of SEQ
    insert_size:    i32,  // INSERT_SIZE tag value, or read_len if insert_size not available
    has_passed_jxn: bool, // READ_HAS_PASSED_JXN tag value
    /* ------------------------------------------- */
    // concatenated alignment-level metadata
    // use String since SvReadPath only exists to be written to file
    node5_observeds:   String, // observed 5' end nodes of the alignment segments
    qry_start0s:       String, // query start positions
    qry_end1s:         String, // query end   positions
    mapqs:             String, // MAPQ values
    cigars:            String, // CIGAR strings
    divergences:       String, // DIVERGENCE tag values
    block_ns:          String, // BLOCK_N tag values
    aln_failure_flags: String, // ALN_FAILURE_FLAG tag values
    jxn_failure_flags: String, // JXN_FAILURE_FLAG tag values (0 or last alignment where not applicable)
    /* ------------------------------------------- */
    // full read SEQ and QUAL as reported by the first alignment of the read
    // reverse complement of the read if seq_strand == -1
    seq_strand: i8, // strand of seq as +1=top|-1=bottom; compare node5_observeds[i].sign() to seq_strand to match seq to cigars[i]
    seq:        String, 
    qual:       String,
}
impl SvReadPath {

    /// Create a new pre-allocated SvReadPath.
    pub fn from_alns(
        alns:    &[BamRecord],
        header:  &HeaderView,
        chroms:  &Chroms,
    ) -> Self{
        let n_alns = alns.len();
        let read_len = alns[0].seq_len();

        // initialize the SvReadPath with read-level values
        let mut read_path = unsafe { SvReadPath{
            qname:            from_utf8_unchecked(alns[0].qname()).to_string(),
            read_len:         read_len as u32,
            insert_size:      tags::get_tag_i32_default( &alns[0], INSERT_SIZE, read_len as i32),
            has_passed_jxn:   tags::get_tag_bool_default(&alns[0], READ_HAS_PASSED_JXN, false),
            /* ------------------------------------------- */
            node5_observeds:  String::with_capacity(n_alns * 14), // try to avoid reallocations
            qry_start0s:      String::with_capacity(n_alns * 10),
            qry_end1s:        String::with_capacity(n_alns * 10),
            mapqs:            String::with_capacity(n_alns * 3),
            cigars:           String::with_capacity(n_alns * 500),
            divergences:      String::with_capacity(n_alns * 20),
            block_ns:         String::with_capacity(n_alns * 3),
            aln_failure_flags:String::with_capacity(n_alns * 4),
            jxn_failure_flags:String::with_capacity(n_alns * 4),
            /* ------------------------------------------- */
            seq_strand:       if alns[0].is_reverse() { -1 } else { 1 },
            seq:              from_utf8_unchecked(&alns[0].seq().as_bytes()).to_string(),
            qual:             from_utf8_unchecked(&alns[0].qual()).to_string(),
        } };

        // iterate to fill in concatenated alignment-level values
        // all values, including single and last values, end with a comma
        alns.iter().for_each(|aln| {
            read_path.node5_observeds  .push_str(&(pack_signed_node_aln(aln, header, 5, chroms).to_string() + COMMA));
            read_path.qry_start0s      .push_str(&(cigar::get_query_start0(aln).to_string() + COMMA));
            read_path.qry_end1s        .push_str(&(cigar::get_query_end1(aln).to_string() + COMMA));
            read_path.mapqs            .push_str(&(aln.mapq().to_string() + COMMA));
            read_path.cigars           .push_str(&(aln.cigar().to_string() + COMMA));
            read_path.divergences      .push_str(&(tags::get_tag_f32(aln, DIVERGENCE).to_string() + COMMA));
            read_path.block_ns         .push_str(&(tags::get_tag_u8_default(aln, BLOCK_N, 1).to_string() + COMMA));
            read_path.aln_failure_flags.push_str(&(tags::get_tag_u8_default(aln, ALN_FAILURE_FLAG, 0).to_string() + COMMA));
            read_path.jxn_failure_flags.push_str(&(tags::get_tag_u8_default(aln, JXN_FAILURE_FLAG, 0).to_string() + COMMA));
        });

        // return the filled SvReadPath
        read_path
    }

    /// Write a vector of AlignmentSegments to a file, one row per
    /// unique sorted segment with a count.
    pub fn write_sorted(
        mut read_paths: Vec<SvReadPath>, 
        bgz_path: &str
    ) -> Result<(), Box<dyn Error>> {
        read_paths.par_sort_unstable();
        let mut bgz_writer = BgzWriter::from_path(bgz_path)?;
        for read_path in &read_paths {
            writeln!(
                bgz_writer, 
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                read_path.qname,
                read_path.read_len,
                read_path.insert_size,
                read_path.has_passed_jxn as u8,
                /* ------------------------------------------- */
                read_path.node5_observeds,
                read_path.qry_start0s,
                read_path.qry_end1s,
                read_path.mapqs,
                read_path.cigars,
                read_path.divergences,
                read_path.block_ns,
                read_path.aln_failure_flags,
                read_path.jxn_failure_flags,
                /* ------------------------------------------- */
                read_path.seq_strand,
                read_path.seq,
                read_path.qual,
            )?;
        }
        Ok(())
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

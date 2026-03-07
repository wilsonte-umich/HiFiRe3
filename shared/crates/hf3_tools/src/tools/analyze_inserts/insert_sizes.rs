//! Report insert sizes and stem lengths for structural variant analysis.

// dependencies
use std::error::Error;
use crossbeam::channel::Sender;
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Counters, Config};
use mdi::OutputFile;
use genomex::sam::SamRecord;
use genomex::genome::Chroms;
use crate::junctions::{JxnFailureFlag};
use crate::formats::hf3_tags::*;
use crate::sites::SiteMatches;
use super::{Outcome, JXN_FAIL_STEM_LENGTH, InsertSizeOutcome, StemLengthOutcome};

// constants
pub_key_constants!{
    // from environment variables
    // IS_COMPOSITE_GENOME
    MIN_SELECTED_SIZE
    SELECTED_SIZE_CV
    MIN_ALLOWED_SIZE
    READ_LENGTH_TYPE
    PLATFORM_MIN_INSERT_SIZE
    PLATFORM_MAX_INSERT_SIZE
    FILTERED_INSERT_SIZES_FILE
    FILTERED_STEM_LENGTHS_FILE
    EXPECTING_ENDPOINT_RE_SITES
    // derived configuration variables
    IS_SIZE_SELECTED
    MAX_ALLOWED_SIZE
    SIZE_PLOT_BIN_SIZE
}
const NON_SV_ACTUAL:     &str = "nonSV.actual";
const NON_SV_PROJECTED:  &str = "nonSV.projected";
const SV_CHIMERIC:       &str = "SV.chimeric"; // cannot assess projected size of structural variant reads
const SV_PASSED_FILTERS: &str = "SV.passed";   // these include all SV reads, not just translocations
// -------------
const INTRAGENOMIC: &str = "intra"; // for assembling types of translocation, i.e., presumed artifacts
const INTERGENOMIC: &str = "inter";
const CHIMERIC:     &str = "chimeric";
const NOT_CHIMERIC: &str = "passed";
const CHIMERIC_DELIMITER: &str = ".";
const SIZE_PLOT_BIN_SIZE_SR: f64 = 10.0;  // short read bin size
const SIZE_PLOT_BIN_SIZE_LR: f64 = 250.0; // long read bin size

/// InsertSizes object for passing insert size outcomes.
pub struct InsertSizes {
    pub insert_size:        i32,
    pub size_bin:           u32,
    pub use_for_size_dist:  bool,
    pub stem_lengths:       Vec<StemLengths>, // one pair of stem lengths per alignment
    pub size_plot_bin_size: f64,
}

/// StemLengths structure for holding 5' and 3' stem lengths.
pub struct StemLengths {
    pub stem5: i32,
    pub stem3: i32,
}

/// InsertSizer helps analyze insert size and junction stem length distributions.
pub struct InsertSizer {
    pub is_size_selected:        bool,
    expecting_endpoint_re_sites: bool,
    min_allowed_size:            i32, // i.e., 1N
    max_allowed_size:            i32, // i.e., 2N
    pub size_plot_bin_size:      f64,
    min_dist_size_bin:           u32,
    max_dist_size_bin:           u32,
}
impl InsertSizer {
    /* ---------------------------------------------------------------------------
    initialize
    ---------------------------------------------------------------------------- */
    /// Initialize an insert size recorder.
    pub fn new(w: &mut Workflow) -> (InsertSizer, Counters, Counters) {
        w.cfg.set_bool_env(  &[EXPECTING_ENDPOINT_RE_SITES]);
        w.cfg.set_string_env(&[READ_LENGTH_TYPE]);
        w.cfg.set_u32_env( &[MIN_SELECTED_SIZE, MIN_ALLOWED_SIZE, PLATFORM_MAX_INSERT_SIZE, PLATFORM_MIN_INSERT_SIZE]);
        w.cfg.set_f64_env(   &[SELECTED_SIZE_CV]);

        let min_allowed_size  = *w.cfg.get_u32(MIN_ALLOWED_SIZE);
        let min_selected_size = *w.cfg.get_u32(MIN_SELECTED_SIZE);
        let selected_size_cv  = *w.cfg.get_f64(SELECTED_SIZE_CV);

        w.cfg.set_bool(IS_SIZE_SELECTED, min_allowed_size > 0 || min_selected_size > 0);
        w.cfg.set_f64(SIZE_PLOT_BIN_SIZE, if w.cfg.equals_string(READ_LENGTH_TYPE, "short") { 
            SIZE_PLOT_BIN_SIZE_SR 
        } else { 
            SIZE_PLOT_BIN_SIZE_LR 
        });
        let is_size_selected   = *w.cfg.get_bool(IS_SIZE_SELECTED);
        let size_plot_bin_size  = *w.cfg.get_f64(SIZE_PLOT_BIN_SIZE);

        if min_allowed_size == 0 {
            w.cfg.set_u32(MIN_ALLOWED_SIZE, (min_selected_size as f64 / (1.0 + selected_size_cv)) as u32);
        }
        let min_allowed_size  = *w.cfg.get_u32(MIN_ALLOWED_SIZE);
        w.cfg.set_u32(MAX_ALLOWED_SIZE, min_allowed_size * 2);
        let max_allowed_size  = *w.cfg.get_u32(MAX_ALLOWED_SIZE);
        let platform_max_insert_size = if is_size_selected {
            max_allowed_size * 4
        } else {
            *w.cfg.get_u32(PLATFORM_MAX_INSERT_SIZE)
        };

        let mut insert_size_counts = Counters::new("insert_size_counts", &[]);
        let mut stem_length_counts = Counters::new("stem_length_counts", &[]);
        insert_size_counts.add_indexed_counters(&[
            (NON_SV_ACTUAL,     "", 0, ""),
            (NON_SV_PROJECTED,  "", 0, ""),
            (SV_CHIMERIC,       "", 0, ""),
            (SV_PASSED_FILTERS, "", 0, ""),
        ]);
        stem_length_counts.add_indexed_counters(&[
            (NON_SV_ACTUAL,     "", 0, ""),
            (NON_SV_PROJECTED,  "", 0, ""),
            (SV_CHIMERIC,       "", 0, ""),
            (SV_PASSED_FILTERS, "", 0, ""),
        ]);

        for translocation_type in [
            format!("{}{}{}", INTRAGENOMIC, CHIMERIC_DELIMITER, CHIMERIC),
            format!("{}{}{}", INTRAGENOMIC, CHIMERIC_DELIMITER, NOT_CHIMERIC),
            format!("{}{}{}", INTERGENOMIC, CHIMERIC_DELIMITER, CHIMERIC),
            format!("{}{}{}", INTERGENOMIC, CHIMERIC_DELIMITER, NOT_CHIMERIC),
        ] {
            insert_size_counts.add_indexed_counters(&[(&translocation_type, "", 0, "")]);
            stem_length_counts.add_indexed_counters(&[(&translocation_type, "", 0, "")]);
        }

        (
            InsertSizer{
                is_size_selected:            is_size_selected,
                expecting_endpoint_re_sites: *w.cfg.get_bool(EXPECTING_ENDPOINT_RE_SITES),
                min_allowed_size:            min_allowed_size as i32,
                max_allowed_size:            max_allowed_size as i32,
                size_plot_bin_size:          size_plot_bin_size,
                min_dist_size_bin:           (*w.cfg.get_u32(PLATFORM_MIN_INSERT_SIZE) as f64 / size_plot_bin_size) as u32,
                max_dist_size_bin:           (platform_max_insert_size                      as f64 / size_plot_bin_size) as u32,
            },
            insert_size_counts,
            stem_length_counts,
        )
        }
    /* ---------------------------------------------------------------------------
    library insert size assessment relative to 1N to 2N expectations
    ---------------------------------------------------------------------------- */
    /// Assess and record insert size and junction stem length for a read pair.
    pub fn assess_insert_sizes(
        &self,
        alns:               &mut [SamRecord],
        is_unmerged_pair:   bool,
        paired_outer_node:  Option<isize>,
        trims:              &[usize],
        is_end_to_end_read: bool,
        chroms:             &Chroms,
    ) -> InsertSizes {
        let n_alns = alns.len();

        // assess insert size, a read-level property
        let mut insert_size = if is_unmerged_pair {
            // when read carries any type of SV, whether sequenced or not, we cannot infer the true insert size
            if n_alns == 1 { // TODO: also check that other alignment does not have a jxn?
                let (_chrom_1, chrom_index_1, ref_pos1_1, is_reverse_1) = 
                    alns[0].get_unpacked_node_aln(5, chroms);
                let (_chrom_2, chrom_index_2, ref_pos1_2, is_reverse_2) = 
                    SamRecord::unpack_signed_node(paired_outer_node.unwrap(), chroms);
                if chrom_index_1 == chrom_index_2 &&
                   is_reverse_1 != is_reverse_2 { // paired outer nodes must have opposite orientations
                    let insert_size = (ref_pos1_2 as i32 - ref_pos1_1 as i32).abs();
                    if insert_size <= self.max_allowed_size * 4 {
                        insert_size
                    } else {
                        self.is_size_selected as i32
                    }
                } else {
                    self.is_size_selected as i32
                }
            } else {
                self.is_size_selected as i32
            }
        } else if alns[0].seq.len() > 1 {
            (alns[0].seq.len() - trims[0] - trims[1]) as i32 // path taken by SV or other reads with sequence
        } else {
            alns[0].get_aligned_size() as i32 // path taken by non-SV reads where SEQ is now *
        };
        let passed_1n_2n = !self.is_size_selected || (
            insert_size >= self.min_allowed_size && 
            insert_size <= self.max_allowed_size
        );
        if !passed_1n_2n { insert_size = -insert_size; }

        // assess junction stem lengths, alignment-level properties
        let qlen = alns[0].seq.len() as u32;
        let query_start0 = alns[0].get_query_start0() as i32;
        let query_end1   = alns[n_alns - 1].get_query_end1(qlen) as i32; // may not be the end of the insert
        let stem_lengths: Vec<StemLengths> = alns.iter().map(|aln| {
            let mut stem5 = aln.get_query_end1(qlen) as i32 - query_start0;
            let passed_stem5 = !self.is_size_selected || stem5 <= self.min_allowed_size;
            if !passed_stem5 { stem5 = -stem5; }
            let stem3 = if is_end_to_end_read {
                let mut stem3 = query_end1 - aln.get_query_start0() as i32;
                let passed_stem3 = !self.is_size_selected || stem3 <= self.min_allowed_size;
                if !passed_stem3 { stem3 = -stem3; }
                stem3
            } else {
                -1 // junctions always fail stem3 test if not end-to-end, where query_end1 != end of insert
            };
            StemLengths { stem5, stem3 }
        }).collect();

        // return results
        let size_bin = (insert_size.abs() as f64 / self.size_plot_bin_size) as u32;
        let use_for_size_dist = size_bin >= self.min_dist_size_bin && size_bin <= self.max_dist_size_bin;
        InsertSizes { 
            insert_size, 
            size_bin, 
            use_for_size_dist, 
            stem_lengths,
            size_plot_bin_size: self.size_plot_bin_size,
        }
    }
    
    /* ---------------------------------------------------------------------------
    perform chimeric junction assessment based on stem length, i.e., <1N logic
    ---------------------------------------------------------------------------- */
    /// Assess whether a junction has fails the stem length filter.
    pub fn check_jxn(
        &self,  
        tx_outcome:   &Sender<Outcome>,
        insert_sizes: &InsertSizes,
        aln5_i:       usize, // the alignment 5' to the junction on the read (its 3' end is the breakpoint)
        aln3_i:       usize,
    ) -> JxnFailureFlag {
        if !self.is_size_selected { return JxnFailureFlag::None; }
        // either breakpoint node can pass the junction stem length test
        let passed = insert_sizes.stem_lengths[aln5_i].stem5 > 0 ||
                           insert_sizes.stem_lengths[aln3_i].stem3 > 0;
        if !passed {
            tx_outcome.send(Outcome::Junction(JXN_FAIL_STEM_LENGTH)).unwrap();
            JxnFailureFlag::StemLength
        } else {
            JxnFailureFlag::None
        }
    }
    /* ---------------------------------------------------------------------------
    export likely SV artifact reads for exploring artifact mechanisms
    ---------------------------------------------------------------------------- */
    /// Record non-SV actual (and projected)sizes.
    pub fn record_non_sv_sizes(
        &self, 
        tx_outcome:   &Sender<Outcome>,
        insert_sizes: &InsertSizes,
        aln_sites:    &[SiteMatches],
    ) -> Result<(), Box<dyn Error>> {
        if !insert_sizes.use_for_size_dist { return Ok(()); }
        tx_outcome.send(Outcome::InsertSize(InsertSizeOutcome{
            read_type:     NON_SV_ACTUAL,
            read_type_str: None,
            size_bin:      insert_sizes.size_bin as usize,
        }))?;
        if !self.expecting_endpoint_re_sites { return Ok(());  }
        let projected_size = (
            aln_sites[0].proj3.pos1 as isize - 
            aln_sites[0].site5.pos1 as isize
        ).abs() + 1;
        let projected_size_bin = (projected_size as f64 / insert_sizes.size_plot_bin_size) as usize;
        tx_outcome.send(Outcome::InsertSize(InsertSizeOutcome{
            read_type:     NON_SV_PROJECTED,
            read_type_str: None,
            size_bin:      projected_size_bin,
        }))?;
        return Ok(()); 
    }
    /// Record likely SV artifact reads for printing to a separate file.
    /// Only on-target reads with exactly two alignments reach this sub.
    pub fn record_sv_sizes(
        &self, 
        tx_outcome:       &Sender<Outcome>,
        insert_sizes:     &InsertSizes,
        alns:             &[SamRecord],
        jxn_failure_flag: u8,
        chroms:           &Chroms,
    ) -> Result<(), Box<dyn Error>> {

        // only tally single-junction SV reads
        // do not record traversal failures or foldback inversions as intermolecular chimeras
        let jxn_flag_init = alns[0].get_tag_value_parsed::<u8>(JXN_FAILURE_FLAG_INIT).unwrap_or(0);
        if alns.len() != 2  || 
           !insert_sizes.use_for_size_dist ||
           JxnFailureFlag::is_pre_chimeric_failure_u8(jxn_flag_init) { 
            return Ok(()); 
        }

        // find the best stem length for this read pair as the shortest distance from
        // the junction to either end of the read (when both are known, o/w use stem5)
        let stem_length5 = insert_sizes.stem_lengths[0].stem5;
        let stem_length3 = insert_sizes.stem_lengths[1].stem3;
        let stem_length3_wrk = if stem_length3.abs() > 1 { stem_length3 } else { -1e9 as i32 }; // only consider meaningful 3' stem lengths
        let min_stem_length = stem_length5.abs().min(stem_length3_wrk.abs());
        let stem_length_bin = (min_stem_length as f64 / self.size_plot_bin_size) as usize;

        // record insert sizes for all single-junction SV reads stratified by chimericity
        // non-chimeric here may contain deletion and other true SVs
        let is_chimeric = JxnFailureFlag::is_chimeric_jxn_u8(jxn_failure_flag);
        if is_chimeric {
            tx_outcome.send(Outcome::InsertSize(InsertSizeOutcome{
                read_type:     SV_CHIMERIC,
                read_type_str: None,
                size_bin: insert_sizes.size_bin as usize,
            }))?;
            tx_outcome.send(Outcome::StemLength(StemLengthOutcome{
                read_type:     SV_CHIMERIC,
                read_type_str: None,
                size_bin:      stem_length_bin,
            }))?;
        } else {
            tx_outcome.send(Outcome::InsertSize(InsertSizeOutcome{
                read_type:     SV_PASSED_FILTERS,
                read_type_str: None,
                size_bin:      insert_sizes.size_bin as usize,
            }))?;
            tx_outcome.send(Outcome::StemLength(StemLengthOutcome{
                read_type:     SV_PASSED_FILTERS,
                read_type_str: None,
                size_bin:      stem_length_bin,
            }))?;
        }

        // only print single-junction interchromosomal chimeras as presumptive artifacts
        // we do not yet know if these are single-molecule or would be confirmed downstream
        // but most interchromosomal, i.e., translocation molecules are likely artifacts
        if alns[0].rname == alns[1].rname { return Ok(()); }

        // stratify interchromosomal molecules as inter- or intra-genomic, when applicable
        let genomicity = if chroms.is_same_genome_suffix(&alns[0].rname, &alns[1].rname) {
            INTRAGENOMIC
        } else {
            INTERGENOMIC
        };

        // reads were stratified above as chimeric or not
        // "chimeric" means a single junction that failed breakpointMatchesSite, isLowQualIns, or junctionHasAdapters
        // "not chimeric" means it passed those chimeric tests, but it may still be chimeric by some other mechanism
        let chimericity = if is_chimeric { CHIMERIC } else { NOT_CHIMERIC };
        let chimera_type = format!("{}{}{}", genomicity, CHIMERIC_DELIMITER, chimericity);

        // increment size counts for translocation presumptive artifacts
        tx_outcome.send(Outcome::InsertSize(InsertSizeOutcome{
            read_type:     "",
            read_type_str: Some(chimera_type.clone()),
            size_bin:      insert_sizes.size_bin as usize,
        }))?;
        tx_outcome.send(Outcome::StemLength(StemLengthOutcome{
            read_type:     "",
            read_type_str: Some(chimera_type),
            size_bin:      stem_length_bin,
        }))?;
        return Ok(()); 
    }

    /// Print insert size and stem length distributions to files for plotting.
    pub fn print_insert_sizes(
        cfg: &mut Config, 
        insert_size_counts: Counters, 
        stem_length_counts: Counters,
        size_plot_bin_size: f64,
    ){
        let mut writer = OutputFile::open_env(cfg, FILTERED_INSERT_SIZES_FILE);
        for key in &insert_size_counts.indexed_keys {
            insert_size_counts.indexed_counts[key].iter().enumerate().for_each(|(bin, count)| {
                let bin = (bin as f64 * size_plot_bin_size) as usize;
                writer.write_record(vec![key, &bin.to_string(), &count.to_string()]);
            });
        }
        writer.close();
        let mut writer = OutputFile::open_env(cfg, FILTERED_STEM_LENGTHS_FILE);
        for key in &stem_length_counts.indexed_keys {
            stem_length_counts.indexed_counts[key].iter().enumerate().for_each(|(bin, count)| {
                let bin = (bin as f64 * size_plot_bin_size) as usize;
                writer.write_record(vec![key, &bin.to_string(), &count.to_string()]);
            });
        }
        writer.close();
    }
}

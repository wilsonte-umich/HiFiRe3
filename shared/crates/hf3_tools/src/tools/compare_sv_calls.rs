//! Compare SVs as called by:
//!     - HiFiRe3 (FinalJunctions) as the test set
//!     - Sniffles2 (VCF) at the truth set
//! Options from `library-properties`, `restriction-enzyme`, and `targets` are
//! applied to the HiFiRe3 test set to establish operating parameters, i.e., to
//! restrict attention to SV junctions from the truth set that could have been 
//! called by the test set. The truth set is used to establish the external
//! list of expected SVs and their RE fragment behavior. Thus, `compare_sv_calls` 
//! supports comparisons of:
//!     - HiFiRe3 and Sniffles2 run on the same data, as well as
//!     - HiFiRe3 run on HiFiRe3 data and Sniffles2 run on external truth data.
//! to benchmark both the HiFiRe3 pipeline alone and the library+pipeline 
//! combination.

// dependencies
use std::error::Error;
use std::fs;
use std::path::PathBuf;
use std::str::from_utf8_unchecked;
// use rustc_hash::FxHashMap;
use rust_htslib::bcf;
use rust_htslib::bcf::Read as VcfRead;
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use mdi::{InputCsv, OutputCsv};
use genomex::genome::{Chroms, TargetRegions};
use crate::junctions::FinalJunction;
use crate::benchmark::*;

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "digest_sv_alleles";
pub_key_constants!(
    // from environment variables
    N_CPU
    SEQUENCING_PLATFORM
    GENOME_FASTA
    IS_COMPOSITE_GENOME
    ENZYME_NAME
    MIN_N_OBSERVED
    MIN_SV_MAPQ
    MIN_SV_LEN
    TRUTHSET_DIR
    TRUTHSET_IS_COMPOSITE_GENOME
    SV_FINAL_JUNCTIONS_FILE_1
    ADMIXED_CALLS_CSV
    // counter keys
    N_INPUT_SVS_HF3
    N_INPUT_SVS_SNF
);

/// Tool for processing VCF SV calls.
struct ComparisonTool {
    // n_cpu: usize,
    min_n_observed: u16,
    min_sv_mapq: u8,
    min_sv_len: u32,
    chroms:  Chroms,
    targets: TargetRegions,
    is_composite_genome: bool,
    truthset_is_composite_genome: bool,
    re: ReSite,
    jxn_calls: Vec<JunctionCall>,
}

// function called by hf3_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[
        SEQUENCING_PLATFORM, GENOME_FASTA, ENZYME_NAME, 
        SV_FINAL_JUNCTIONS_FILE_1, TRUTHSET_DIR,
        ADMIXED_CALLS_CSV
    ]);
    cfg.set_usize_env( &[N_CPU, MIN_N_OBSERVED, MIN_SV_MAPQ, MIN_SV_LEN]);
    cfg.set_bool_env( &[IS_COMPOSITE_GENOME, TRUTHSET_IS_COMPOSITE_GENOME]);

    // initialize counters
    let ctrs = Counters::new(TOOL, &[
        (N_INPUT_SVS_HF3,     "number of input SV lines in HiFiRe3 final junctions"),
        (N_INPUT_SVS_SNF,     "number of input SV lines in truth-set VCF"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();

    // collect the working chromosomes and targets
    let chroms = Chroms::new(&mut w.cfg);
    let targets = TargetRegions::from_env(&mut w, false);

    // create the digestion tool
    let mut tool = ComparisonTool {
        // n_cpu:          *w.cfg.get_usize(N_CPU),
        min_n_observed: *w.cfg.get_usize(MIN_N_OBSERVED) as u16,
        min_sv_mapq:    *w.cfg.get_usize(MIN_SV_MAPQ) as u8,
        min_sv_len:     *w.cfg.get_usize(MIN_SV_LEN) as u32,
        chroms:  chroms,
        targets: targets,
        is_composite_genome: *w.cfg.get_bool(IS_COMPOSITE_GENOME),
        truthset_is_composite_genome: *w.cfg.get_bool(TRUTHSET_IS_COMPOSITE_GENOME),
        re: ReSite::from_re_name(w.cfg.get_string(ENZYME_NAME), None),
        jxn_calls: Vec::new(),
    };

    // load junctions into a single table identified by source
    load_hf3_jxns(&mut w, &mut tool)?;
    load_snf_jxns(&mut w, &mut tool)?;
    write_admixed_calls(&mut w, &mut tool);

    // print counts
    w.ctrs.print_grouped(&[
        &[N_INPUT_SVS_HF3, N_INPUT_SVS_SNF],
    ]);
    Ok(())
}

// -----------------------------------------------------------------------------
// load junctions from both data sets
// -----------------------------------------------------------------------------
/// Load HiFiRe3 Final Junctions.
fn load_hf3_jxns(w: &mut Workflow, tool: &mut ComparisonTool) -> Result<(), Box<dyn Error>>{
    w.log.print("loading HiFiRe3 junctions");
    let filepath = w.cfg.get_string(SV_FINAL_JUNCTIONS_FILE_1);
    let mut csv = InputCsv::open(filepath);
    for result in csv.deserialize::<FinalJunction>() {
        let jxn = result?;
        let mut qnames: Vec<_> = jxn.qnames.split(",").collect();
        qnames.sort_unstable();
        let qnames = qnames.join(",");

        // only continue with primary species canonical chromosomes
        let Some(chrom_1) = hf3_chrom_index_to_name(
            jxn.chrom_index1_1, tool.is_composite_genome, &tool.chroms
        ) else { continue };
        let Some(chrom_2) = hf3_chrom_index_to_name(
            jxn.chrom_index1_2, tool.is_composite_genome, &tool.chroms
        ) else { continue };

        // collect additional values
        let n_observed = jxn.n_instances_dedup;

        // apply additional filters
        if n_observed < tool.min_n_observed { continue }
        if n_observed < tool.min_n_observed { continue }
        if n_observed < tool.min_n_observed { continue }
        if n_observed < tool.min_n_observed { continue }
        if n_observed < tool.min_n_observed { continue }

        tool.jxn_calls.push(JunctionCall { 
            source:     SvSource::HiFiRe3,
            sv_type:    SvType::from_final_junction(&jxn), 
            sv_size:    jxn.sv_size, 
            chrom_1, 
            pos1_1:     jxn.ref_pos1_1, 
            chrom_2, 
            pos1_2:     jxn.ref_pos1_2, 
            jxn_seq:    jxn.jxn_seq.clone(),
            n_observed,
            frag_size:  jxn.get_median_insert_size() as i32,
            qnames,
        });
    }
    Ok(())
}

// shared:
//         min-selected-size: 0
//         min-sv-mapq:    50

fn load_snf_jxns(w: &mut Workflow, tool: &mut ComparisonTool) -> Result<(), Box<dyn Error>>{
    w.log.print("loading truth-set junctions");

    // find the truth set file
    let dir = w.cfg.get_string(TRUTHSET_DIR);
    let mut filepath = None::<PathBuf>;
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_file() {
            if path.ends_with("frag_sizes.vcf.gz"){
                filepath = Some(path);
            }
        }
    }
    let filepath = filepath.expect("missing truth set VCF");

    // process SV lines (which may match two inversion junctions)
    let mut vcf = bcf::Reader::from_path(filepath)?;
    for result in vcf.records() {
        let sv = result?;

        // only continue with primary species canonical chromosomes
        let alt  = get_vcf_alt(&sv)?;
        let (chrom_1, pos1_1) = get_vcf_left_node(&sv)?;
        let (chrom_2, pos1_2) = get_vcf_right_node(&sv, &alt, &chrom_1)?;
        let Some(chrom_1) = simplify_chrom_name(
            &chrom_1, tool.truthset_is_composite_genome
        ) else { continue };
        let Some(chrom_2) = simplify_chrom_name(
            &chrom_2, tool.truthset_is_composite_genome
        ) else { continue };
        let pos1_1 = pos1_1 as u32;
        let pos1_2 = pos1_2 as u32;

        // filter to retain only on-target junctions
        let is_on_target = !tool.targets.has_data || 
            tool.targets.get_pos_region(&chrom_1, pos1_1).2 == INSIDE_REGION ||
            tool.targets.get_pos_region(&chrom_2, pos1_2).2 == INSIDE_REGION;
        if !is_on_target { continue }

        // collect SV-level values, the same for all junctions
        let sv_type = SvType::from_vcf(&sv);
        let sv_size = get_vcf_sv_len_with_default(&sv)?;
        let n_observed = get_vcf_n_supplemental(&sv)?;
        let frag_size = match sv.info(&tool.re.size_tag).integer()? {
            Some(v) => (v[0] as i32).to_owned(),
            None => 0 // the path taken when RE is NA
        };

        // apply additional filters
        if n_observed < tool.min_n_observed { continue }
        if sv_size < tool.min_sv_len { continue }

        let qnames = get_vcf_qnames(&sv)?.expect("missing qnames");
        let mut qnames: Vec<_> = qnames.iter().map(|x|{
            unsafe{ from_utf8_unchecked(x) }
        }).collect();
        qnames.sort_unstable();
        let qnames = qnames.join(",");

        match sv_type {
            SvType::Inv | SvType::Ins => {
                // TODO: handle later
            },
            _ => tool.jxn_calls.push(JunctionCall { 
                source: SvSource::TruthSet,
                sv_type, 
                sv_size, 
                chrom_1, 
                pos1_1, 
                chrom_2, 
                pos1_2, 
                jxn_seq: get_vcf_ins_seqs(&sv, 0)?,
                n_observed,
                frag_size,
                qnames
            })
        };
    }
    Ok(())
}

// -----------------------------------------------------------------------------
// write a table of admixed SV/junction calls
// -----------------------------------------------------------------------------
fn write_admixed_calls(w: &mut Workflow, tool: &mut ComparisonTool){
    w.log.print("writing admixed calls");
    let csv = OutputCsv::open_env(&mut w.cfg, ADMIXED_CALLS_CSV, None);
    csv.serialize_all::<JunctionCall>(&tool.jxn_calls);
}

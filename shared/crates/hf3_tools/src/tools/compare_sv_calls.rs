//! Compare SVs as called by:
//!     - HiFiRe3 (FinalJunctions) as the test set
//!     - Sniffles2 and CuteSV (VCF) at two truth sets
//! Options from `library-properties`, `restriction-enzyme`, and `targets` are
//! applied to the HiFiRe3 test set to establish operating parameters, i.e., to
//! restrict attention to SV junctions from the truth sets that could have been 
//! called by the test set. The truth sets are used to establish the external
//! list of expected SVs and their RE fragment behavior. Thus, `compare_sv_calls` 
//! supports comparisons of:
//!     - HiFiRe3 and Sniffles2/CuteSV run on the same data, as well as
//!     - HiFiRe3 run on HiFiRe3 data and Sniffles2/CuteSV run on external truth data.
//! to benchmark both the HiFiRe3 pipeline alone and the library+pipeline 
//! combination. Along the way, Sniffles2 and CuteSV are also compared.

// dependencies
use std::error::Error;
use std::fs;
use std::path::PathBuf;
// use std::str::from_utf8_unchecked;
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
const TOOL: &str = "compare_sv_calls";
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
pub struct ComparisonTool {
    // n_cpu: usize,
    min_n_observed: u16,
    _min_sv_mapq: u8,
    min_sv_len: u32,
    chroms:  Chroms,
    targets: TargetRegions,
    is_composite_genome: bool,
    truthset_is_composite_genome: bool,
    pub re: ReSite,
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
        _min_sv_mapq:    *w.cfg.get_usize(MIN_SV_MAPQ) as u8,
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
    load_vcf_jxns(&mut w, &mut tool, "sniffles")?;
    load_vcf_jxns(&mut w, &mut tool, "cutesv")?;
    write_admixed_calls(&mut w, &mut tool);

    // print counts
    w.ctrs.print_grouped(&[
        &[N_INPUT_SVS_HF3, N_INPUT_SVS_SNF],
    ]);
    Ok(())
}

// -----------------------------------------------------------------------------
// load junctions from all data sets
// -----------------------------------------------------------------------------
/// Load HiFiRe3 Final Junctions.
fn load_hf3_jxns(w: &mut Workflow, tool: &mut ComparisonTool) -> Result<(), Box<dyn Error>>{
    w.log.print("loading HiFiRe3 junctions");

    // open the HiFiRe3 final junctions file
    let filepath = w.cfg.get_string(SV_FINAL_JUNCTIONS_FILE_1);
    let mut csv = InputCsv::open(filepath);

    // process junction lines (complete inversion have two associated junction lines)
    for result in csv.deserialize::<FinalJunction>() {
        let jxn = result?;

        // only continue with primary species canonical chromosomes
        let Some(chrom_1) = hf3_chrom_index_to_name(
            jxn.chrom_index1_1, tool.is_composite_genome, &tool.chroms
        ) else { continue };
        let Some(chrom_2) = hf3_chrom_index_to_name(
            jxn.chrom_index1_2, tool.is_composite_genome, &tool.chroms
        ) else { continue };

        // filter to retain only on-target junctions
        // redundant with SV calling, but ensure consistency with load_snf_jxns()
        let mut pos1_1 = jxn.ref_pos1_1;
        let mut pos1_2 = jxn.ref_pos1_2; 
        if !sv_is_on_target(&chrom_1, pos1_1, &chrom_1, pos1_1, &tool.targets) { 
            continue 
        }

        // collect additional values
        let n_observed = jxn.n_instances_dedup; //n_reads_dedup

        // apply additional filters
        if n_observed < tool.min_n_observed { continue }
        if jxn.sv_size < tool.min_sv_len { continue }

        // adjust output to be as consistent with VCF encoding
        let sv_type = SvType::from_final_junction(&jxn);
        let mut sv_size = jxn.sv_size;
        match sv_type {
            SvType::Del => {

// chr1_hs1        233212921       Sniffles2.DEL.7FS0      g       <DEL>   60      PASS    PRECISE;SVTYPE=DEL;SVLEN=-1341;END=233214262;SUPPORT=6;RNAMES=m21026_251212_225317/145884959/ccs/rev,m21026_251212_225317/152176957/ccs/fwd,m21026_251212_225317/186582474/ccs/rev,m21026_251212_225317/145884959/ccs/fwd,m21026_251212_225317/152176957/ccs/rev,m21026_251212_225317/186582474/ccs/fwd;COVERAGE=0,0,0,0,6;STRAND=+-;NM=0.008;PHASE=0,NULL,6,6,PASS,FAIL;STDEV_LEN=0.000;STDEV_POS=0.000;SUPPORT_SA=6;VAF=1.000       GT:GQ:DR:DV:PS  1/1:16:0:6:.

// chr1_hs1        233212921       cuteSV.DEL.287  gt...gt g       .       PASS    PRECISE;SVTYPE=DEL;SVLEN=-1342;END=233214263;CIPOS=-0,0;CILEN=-0,0;RE=6;RNAMES=m21026_251212_225317/145884959/ccs/rev,m21026_251212_225317/145884959/ccs/fwd,m21026_251212_225317/152176957/ccs/fwd,m21026_251212_225317/152176957/ccs/rev,m21026_251212_225317/186582474/ccs/fwd,m21026_251212_225317/186582474/ccs/rev;STRAND=+-  GT:DR:DV:PL:GQ  ./.:.:6:.,.,.:.

// del     chr1    233212921       chr1    233214262       **      SNF     1341    3743    513     0       TT
// del     chr1    233212921       chr1    233214263       +-      CSV     1342    3743    514     0       TG
// del     chr1    233212921       chr1    233214264       +-      HF3     1343    3745    516     -2      *

                // -----MMM-----
                // -----MMM~~~~~~~MMM-----
                // H      1       2
                // V      1         2
                // -----|-----
                // -----~~~~~~~-----
                // H   1       2
                // V   1      2
                if jxn.offset < 0 { pos1_2 += jxn.offset.unsigned_abs() as u32 }
                pos1_2 -= 1;
                sv_size = pos1_2 - pos1_1;
            },
            SvType::Dup => {
                // -----MMM~~~~~~MMM~~~~~~MMM-----
                // -----MMM~~~~~~MMM-----
                // H    2          1
                // V      1        2
                // -----~~~~~~~|~~~~~~~-----
                // -----~~~~~~~-----
                // H    2     1
                // V   1      2
                if jxn.offset < 0 { pos1_2 += jxn.offset.unsigned_abs() as u32 }
                pos1_2 -= 1;
                (pos1_1, pos1_2) = (pos1_2, pos1_1);
                sv_size = pos1_2 - pos1_1;
            },
            SvType::Inv if jxn.strand_index0_1 == 0 => {
                // -----MMM7654321MMM-----
                // -----MMM1234567MMM-----
                // 
                // -----|7654321|-----
                // -----1234567-----
                // H   1      2
                // V   1  2
            },
            SvType::Inv if jxn.strand_index0_1 == 1 => {
                // -----|7654321|-----
                // -----1234567-----
                // H    1      2
                // V   
            },
            _ => {},
        }

        // assemble a sorted qname concatenation
        let qnames = get_sorted_qnames_str(&jxn.qnames);
        // commit the junction call
        tool.jxn_calls.push(JunctionCall { 
            sv_type, 

            chrom_1, 
            pos1_1, 
            chrom_2, 
            pos1_2, 
            strands: jxn.get_strand_encoding(),

            source:     SvSource::HiFiRe3,

            sv_size:     sv_size, 
            frag_size:   jxn.get_median_insert_size() as i32,
            stem_length: jxn.min_stem_length,

            offset:  jxn.offset,
            ins_seq: if jxn.offset > 0 { jxn.jxn_seq.clone() } else { "*".to_string() },
            
            n_observed,
            mapq:       jxn.max_min_mapq,
            divergence: jxn.min_max_divergence,

            jxn_fail:     jxn.jxn_failure_flag,
            aln_fail:     jxn.aln_failure_flag,
            all_chimeric: !jxn.any_was_not_chimeric,

            mosaic: false,
            vaf: 0.0,

            qnames,
        });
    }
    Ok(())
}

/// Load truth-set VCF as generated by benchmark::align+digest.
fn load_vcf_jxns(w: &mut Workflow, tool: &mut ComparisonTool, sv_caller: &str) -> Result<(), Box<dyn Error>>{
    w.log.print(&format!("loading {} truthset junctions", sv_caller));
    let sv_source = match sv_caller {
        "sniffles" => SvSource::Sniffles,
        _ =>          SvSource::CuteSV,
    };

    // find and open the truth set VCF file
    let dir = w.cfg.get_string(TRUTHSET_DIR);
    let mut filepath = None::<PathBuf>;
    for entry in fs::read_dir(dir)? {
        let path = entry?.path();
        if path.is_file() && path
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.ends_with(&format!("{}.frag_sizes.vcf.gz", sv_caller))) {
            filepath = Some(path);
            break;
        }
    }
    let filepath = filepath.expect(&format!("missing {} truth set VCF", sv_caller));
    let mut vcf = bcf::Reader::from_path(filepath)?;

    // process junction lines, as filtered and modified by `digest_sv_alleles`
    for result in vcf.records() {
        let jxn = result?;
        let alt = get_vcf_alt(&jxn)?;
        let (sv_type, sv_type_vcf) = SvType::from_vcf(&jxn);

        // only continue with primary species canonical chromosomes
        let (chrom_1, pos1_1) = get_vcf_left_node(&jxn)?;
        let (chrom_2, pos1_2) = get_vcf_right_node(&jxn, &chrom_1, &alt, &sv_type_vcf)?;
        let Some(chrom_1) = simplify_chrom_name(
            &chrom_1, tool.truthset_is_composite_genome
        ) else { continue };
        let Some(chrom_2) = simplify_chrom_name(
            &chrom_2, tool.truthset_is_composite_genome
        ) else { continue };

        // filter to retain only on-target junctions
        if !sv_is_on_target(&chrom_1, pos1_1, &chrom_1, pos1_1, &tool.targets) { 
            continue 
        }

        // collect additional values
        let sv_size = get_vcf_sv_size(&chrom_1, pos1_1, &chrom_2, pos1_2);
        let n_observed = get_vcf_n_observed(&jxn, &sv_source)?;
        let qnames = get_vcf_qnames(&jxn)?.expect("Missing INFO/RNAMES");

        // apply additional filters
        if n_observed < tool.min_n_observed { continue }
        if sv_size < tool.min_sv_len { continue }

        // commit junction call
        tool.jxn_calls.push(JunctionCall { 
            sv_type, 

            chrom_1, 
            pos1_1, 
            chrom_2, 
            pos1_2, 
            strands: get_vcf_strand_str(&jxn, &sv_source)?,

            source: sv_source.clone(),

            sv_size, 
            frag_size:   get_frag_size(&jxn, &tool)?,
            stem_length: get_stem_length(&jxn, &tool)?,

            offset: 0,
            ins_seq: get_vcf_ins_seq(&jxn)?,
            
            n_observed,
            mapq: jxn.qual() as u8,
            divergence: 0.0,

            jxn_fail: 0,
            aln_fail: 0,
            all_chimeric: 0,

            mosaic: get_vcf_mosaic(&jxn, &sv_source)?,
            vaf:    get_vcf_vaf(&jxn, &sv_source)?,

            qnames: get_sorted_qnames_bytes(qnames),
        });
    }
    Ok(())
}

// -----------------------------------------------------------------------------
// write a table of admixed SV/junction calls
// -----------------------------------------------------------------------------
fn write_admixed_calls(w: &mut Workflow, tool: &mut ComparisonTool){
    w.log.print("sorting and writing admixed calls");
    tool.jxn_calls.sort_unstable();
    let csv = OutputCsv::open_env(&mut w.cfg, ADMIXED_CALLS_CSV, None);
    csv.serialize_all::<JunctionCall>(&tool.jxn_calls);
}

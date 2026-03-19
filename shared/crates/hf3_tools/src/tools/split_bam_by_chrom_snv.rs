//! Split input name-sorted BAM file(s) into temporary per-chromosome BAM 
//! files based on the chromosome of each alignment. Unlike for SVs,
//! alignments from multi-alignment reads may be dispatched to different 
//! files.
//! 
//! Output only includes usable on-target PacBioStrand reads, as defined 
//! by the presence of the minimap2 cs tag.
//! 
//! Support multiple BAM files for multi-sample variant calling.

// dependencies
use std::error::Error;
use std::str::from_utf8_unchecked;
use rustc_hash::FxHashMap;
use rust_htslib::bam::{Reader, Read, Writer, Record as BamRecord, Header, HeaderView, Format, record::Aux};
use rust_htslib::tpool::ThreadPool;
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use mdi::OutputFile;
use genomex::genome::{Chroms, TargetRegions};
use genomex::bam::qual::median_qual_aln;
use crate::formats::hf3_tags::*;
use crate::snvs::check_pacbio_strand;

// constants
const TOOL: &str = "split_by_chrom_snv";
pub_key_constants!(
    // from environment variables
    N_CPU
    IS_COMPOSITE_GENOME
    MIN_AVG_BASE_QUAL
    NAME_BAM_FILES
    INDEX_FILE_PREFIX_WRK
    SNV_SAMPLES_FILE
    // counter keys
    N_ALNS
    N_USABLE_ALNS
    N_BASES
    N_ALNS_BY_GENOME
    N_BASES_BY_GENOME
    N_ALNS_BY_SAMPLE
    N_BASES_BY_SAMPLE
);
const CS_TAG: &[u8] = DIFFERENCE_STRING.as_bytes();
const SB_TAG: &[u8] = SAMPLE_BIT.as_bytes();
const MIN_MAPQ: u8  = 50; // TODO: expose as option if also done so for strand_merger

// main function called by xxx_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_u32_env(&[N_CPU]);
    cfg.set_bool_env(&[IS_COMPOSITE_GENOME]);
    cfg.set_u8_env(&[MIN_AVG_BASE_QUAL]);
    cfg.set_string_env(&[NAME_BAM_FILES, INDEX_FILE_PREFIX_WRK, SNV_SAMPLES_FILE]);
    let min_avg_base_qual = *cfg.get_u8(MIN_AVG_BASE_QUAL);

    // validate we are working with the expected read data type
    check_pacbio_strand(TOOL, &mut cfg)?;

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_ALNS,        "alignments processed"),
        (N_USABLE_ALNS, "usable on-target alignments in output"),
        (N_BASES,       "reference bases in usable on-target alignments"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_ALNS_BY_GENOME,  "on-target alignments by genome"),
        (N_BASES_BY_GENOME, "reference bases in on-target alignments by genome"),
        (N_ALNS_BY_SAMPLE,  "on-target alignments by sample"),
        (N_BASES_BY_SAMPLE, "reference bases in on-target alignments by sample"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();

    // collect the working chromosomes
    let chroms = Chroms::new(&mut w.cfg);
    let targets = TargetRegions::from_env(&mut w, false);
    let on_target_chroms = targets.get_region_chroms(&chroms);
    let is_composite = *w.cfg.get_bool(IS_COMPOSITE_GENOME);
    let chrom_file_prefix = w.cfg.get_string(INDEX_FILE_PREFIX_WRK);

    // use a thread pool for BAM reading and writing
    let tpool = ThreadPool::new(w.cfg.get_u32(N_CPU) - 1)?;

    // initialize the output BAM writers, one per target chromosome shared over all samples
    let name_bam_files = w.cfg.get_string(NAME_BAM_FILES);
    let name_bam_paths = name_bam_files.split(',').collect::<Vec<&str>>();
    let name_bam = Reader::from_path(name_bam_paths[0])?;
    let header_view = name_bam.header().clone(); // for TID lookups
    let header = Header::from_template(&header_view); // shared header for each output BAM writer
    let mut writers:     FxHashMap<u32, Writer> = FxHashMap::default(); // TID -> file writer named by our padded chrom index
    for (chrom, chrom_index) in on_target_chroms.iter() {
        let tid = header_view
            .tid(chrom.as_bytes())
            .expect(format!("{} not found in BAM header", chrom).as_str()) as u32;
        let chrom_index_padded = format!("{:02}", chrom_index);
        let mut writer = Writer::from_path(
            format!("{}.chr{}.bam", chrom_file_prefix, chrom_index_padded),
            &header,
            Format::Bam
        ).expect(&format!("Failed to create BAM writer for chrom {}", chrom));
        writer.set_thread_pool(&tpool)?;
        writers.insert(tid, writer);
    }
    drop(name_bam);

    // run through multiple name BAM files to support single and multi-sample analyses
    w.log.print("streaming BAM records");
    let samples_file = w.cfg.get_string(SNV_SAMPLES_FILE);
    let header = vec!["sample_bit", "sample_name"];
    let mut samples_file = OutputFile::open_file(&samples_file, b'\t', Some(&header)); 
    let mut sample_bit: u32 = 1;
    for name_bam_path in name_bam_paths {
        let bam_file_name = name_bam_path.split('/').last().unwrap();
        let sample_name = bam_file_name.split('.').nth(0).unwrap();
        samples_file.write_record(vec![&sample_bit.to_string(), sample_name]);
        let mut name_bam = Reader::from_path(name_bam_path)?;
        name_bam.set_thread_pool(&tpool)?;

        // process input BAM records
        eprintln!("    {}", sample_name);
        let mut aln = BamRecord::new();
        while let Some(result) = name_bam.read(&mut aln) {
            match result {
                Ok(_)  => print_aln(
                    &mut aln, 
                    &header_view, 
                    is_composite, 
                    min_avg_base_qual,
                    &mut writers, 
                    &mut w.ctrs, 
                    sample_bit, 
                    sample_name
                )?,
                Err(_) => panic!("BAM parsing failed")
            }
        }
        sample_bit <<= 1; 
    }   
    
    // report counter values
    w.ctrs.print_grouped(&[
        &[N_ALNS, N_USABLE_ALNS, N_BASES],
        &[N_ALNS_BY_GENOME],
        &[N_BASES_BY_GENOME],
        &[N_ALNS_BY_SAMPLE],
        &[N_BASES_BY_SAMPLE],
    ]);
    Ok(())
}

// print and count on-target alignments
fn print_aln(
    aln:          &mut BamRecord, 
    header_view:  &HeaderView,
    is_composite: bool,
    min_avg_base_qual: u8,
    writers:      &mut FxHashMap<u32, Writer>,
    ctrs:         &mut Counters,
    sample_bit:   u32,
    sample_name:  &str,
) -> Result<(), Box<dyn Error>> {
    ctrs.increment(N_ALNS);

    // skip unmapped reads
    // rare mapped reads lack a cs tag, so the presence of a cs tag is the most specific criterion
    if let Some(_cs) = aln.aux(CS_TAG).ok() { 

        // skip low quality alignments
        if aln.mapq() < MIN_MAPQ || 
           median_qual_aln(aln) < min_avg_base_qual { 
            return Ok(()); 
        }

        // skip reads in untargeted samples that map to other than nuclear chromosomes
        let tid = aln.tid();
        if let Some(writer) = writers.get_mut(&(tid as u32)){
            ctrs.increment(N_USABLE_ALNS);
            ctrs.increment_keyed(N_ALNS_BY_SAMPLE, sample_name);

            // add the sample bit tag for multi-sample comparison
            aln.push_aux(SB_TAG, Aux::U32(sample_bit)).unwrap();

            // commit on-target reads to temporary BAM files
            writer.write(aln)?; 

            // increment counters
            let cigar_view = aln.cigar();
            let n_bases = cigar_view.end_pos() as usize - aln.pos() as usize;
            ctrs.add_to(N_BASES, n_bases); 
            ctrs.add_to_keyed(N_BASES_BY_SAMPLE, sample_name, n_bases);
            if is_composite {
                let chrom = unsafe{ from_utf8_unchecked(header_view.tid2name(aln.tid() as u32))}; 
                let genome_name = chrom.split_once('_').unwrap().1; // e.g., chr1_hs1
                ctrs.increment_keyed(N_ALNS_BY_GENOME, genome_name);
                ctrs.add_to_keyed(N_BASES_BY_GENOME,  genome_name, n_bases);
            }
        }
    }
    Ok(())
}

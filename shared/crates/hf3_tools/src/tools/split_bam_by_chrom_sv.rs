//! Split input name-sorted BAM file(s) into temporary per-chromosome BAM 
//! files based on the chromosome of the first alignment in each read.
//! 
//! Output only includes usable on-target reads, as defined by the absence 
//! the UNMAPPED flag bit and READ_IS_OFF_TARGET and READ_FAILURE_FLAG tags.
//! 
//! Along the way, collect de-duplicated on-target coverage statistics.
//! 
//! Support multiple BAM files for multi-sample variant calling.

// dependencies
use std::error::Error;
use std::fs::File;
use std::io::Write;
use rustc_hash::FxHashMap;
use rayon::prelude::*;
use rust_htslib::bam::{Reader, Read, Writer, Record as BamRecord, Header, Format, record::Aux};
use rust_htslib::tpool::ThreadPool;
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use mdi::OutputFile;
use genomex::genome::{Chroms, TargetRegions};
use genomex::sam::SamRecord;
use crate::formats::hf3_tags::*;
use crate::inserts::{UniqueInsertSpan, ChannelAlignment};

// constants
const TOOL: &str = "split_by_chrom";
pub_key_constants!(
    // from environment variables
    N_CPU
    DEDUPLICATE_READS
    IS_COMPOSITE_GENOME
    NAME_BAM_FILES
    INDEX_FILE_PREFIX_WRK
    SV_SAMPLES_FILE
    // counter keys
    N_READS
    N_USABLE_READS
    N_SV_READS
    N_UNIQ_ALNS
    N_ALNS
    N_BASES
    N_READS_BY_GENOME
    N_BASES_BY_GENOME
    N_READS_BY_SAMPLE
    N_BASES_BY_SAMPLE
);
const READ_FAILURE: &[u8] = READ_FAILURE_FLAG.as_bytes();
const OFF_TARGET: &[u8]   = READ_IS_OFF_TARGET.as_bytes(); 
const SB_TAG: &[u8] = SAMPLE_BIT.as_bytes();
const COUNTER_CAPACITY: usize = 100_000_000;

// main function called by xxx_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_u32_env(&[N_CPU]);
    cfg.set_bool_env(&[DEDUPLICATE_READS, IS_COMPOSITE_GENOME]);
    cfg.set_string_env(&[NAME_BAM_FILES, INDEX_FILE_PREFIX_WRK, SV_SAMPLES_FILE]);

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_READS,        "reads processed"),
        (N_USABLE_READS, "usable on-target reads in output"),
        (N_SV_READS,     "usable on-target reads with at least one SV junction in output"),
        (N_UNIQ_ALNS,    "unique (deduplicated) alignments in on-target reads"),
        (N_ALNS,         "total  (deduplicated) alignments in on-target reads"),
        (N_BASES,        "(deduplicated) reference bases in on-target alignments"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_READS_BY_GENOME, "on-target reads by genome"),
        (N_BASES_BY_GENOME, "reference bases in on-target alignments by genome"),
        (N_READS_BY_SAMPLE, "on-target reads by sample"),
        (N_BASES_BY_SAMPLE, "reference bases in on-target alignments by sample"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();

    // collect the working chromosomes
    let chroms = Chroms::new(&mut w.cfg);
    let targets = TargetRegions::from_env(&mut w, false);
    let on_target_chroms = targets.get_region_chroms(&chroms);

    // use a thread pool for BAM reading and writing
    let tpool = ThreadPool::new(w.cfg.get_u32(N_CPU) - 1)?;

    // initialize the output BAM writers, one per target chromosome    let name_bam_files = w.cfg.get_string(NAME_BAM_FILES);
    let name_bam_files = w.cfg.get_string(NAME_BAM_FILES);
    let name_bam_paths = name_bam_files.split(',').collect::<Vec<&str>>();
    let name_bam = Reader::from_path(name_bam_paths[0])?;
    let chrom_file_prefix = w.cfg.get_string(INDEX_FILE_PREFIX_WRK);
    let header_view = name_bam.header(); // for TID lookups
    let header = Header::from_template(header_view); // shared header for each output BAM writer
    let mut writers:     FxHashMap<u32, Writer> = FxHashMap::default(); // TID -> file writer named by our padded chrom index
    let mut read_counts: FxHashMap<u32, usize>  = FxHashMap::default(); // TID -> count of reads
    let mut chrom_map:   FxHashMap<u32, String> = FxHashMap::default(); // TID -> padded chrom index
    for (chrom, chrom_index) in on_target_chroms {
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
        read_counts.insert(tid, 0); 
        chrom_map.insert(tid, chrom_index_padded);
    }
    drop(name_bam);

    // initialize the deduplication counter
    let mut channel_alns: Vec<(ChannelAlignment, u32)> = Vec::with_capacity(COUNTER_CAPACITY);

    // run through multiple name BAM files to support single and multi-sample analyses
    w.log.print("streaming BAM records");
    let samples_file = w.cfg.get_string(SV_SAMPLES_FILE);
    let header = vec!["sample_bit", "sample_name"];
    let mut samples_file = OutputFile::open_file(&samples_file, b'\t', Some(&header)); 
    let mut sample_bit: u32 = 1;
    let mut samples: FxHashMap<u32, String> = FxHashMap::default(); // sample_bit -> sample_name
    for name_bam_path in name_bam_paths {
        let bam_file_name = name_bam_path.split('/').last().unwrap();
        let sample_name = bam_file_name.split('.').nth(0).unwrap();
        samples_file.write_record(vec![&sample_bit.to_string(), sample_name]);
        samples.insert(sample_bit, sample_name.to_string());
        let mut name_bam = Reader::from_path(name_bam_path)?;
        name_bam.set_thread_pool(&tpool)?;

        // process input BAM records
        eprintln!("    {}", sample_name);
        let mut records = name_bam.records();
        let mut alns: Vec<BamRecord> = vec![records.next().unwrap()?]; // expect there to always be data, thus alns is never empty
        for result in records {
            let record = result?;
            if record.qname() != alns[0].qname() { 
                print_alns(
                    &mut alns, 
                    &mut writers, &mut read_counts, &mut w.ctrs, &mut channel_alns, 
                    sample_bit, sample_name
                )?;
                alns.clear();
            }
            alns.push(record);
        }
        print_alns(
            &mut alns,
            &mut writers, &mut read_counts, &mut w.ctrs, &mut channel_alns,
            sample_bit, sample_name
        )?;
        sample_bit <<= 1; 
    }

    // write chrom read count files, for setting capacity in downstream tools
    w.log.print("writing chromosome read count files");
    for (tid, count) in &read_counts {
        let chrom_index_padded = chrom_map.get(tid).unwrap();
        let count_file_path = format!("{}.chr{}.read_count", chrom_file_prefix, chrom_index_padded);
        let mut count_file = File::create(&count_file_path)?;
        writeln!(count_file, "{}", count)?;
    }
    
    // report counter values
    w.log.print("calling tally_read_bases");
    tally_read_bases(&mut channel_alns, &chroms, &mut w, samples);
    w.ctrs.print_grouped(&[
        &[N_READS, N_USABLE_READS, N_SV_READS],
        &[N_UNIQ_ALNS, N_ALNS, N_BASES],
        &[N_READS_BY_GENOME],
        &[N_BASES_BY_GENOME],
        &[N_READS_BY_SAMPLE],
        &[N_BASES_BY_SAMPLE],
    ]);
    Ok(())
}

// print and count on-target alignments
fn print_alns(
    alns:         &mut [BamRecord], 
    writers:      &mut FxHashMap<u32, Writer>,
    read_counts:  &mut FxHashMap<u32, usize>,
    ctrs:         &mut Counters,
    channel_alns: &mut Vec<(ChannelAlignment, u32)>,
    sample_bit:   u32,
    sample_name:  &str,
) -> Result<(), Box<dyn Error>> {

    // skip unmapped, off-target, or failed reads
    ctrs.increment(N_READS);
    let tid = alns[0].tid();
    if tid < 0 || // skip unmapped reads (synonymous with UNMAPPED flag bit being set)
       alns[0].aux(READ_FAILURE).is_ok() || // skip if read failure flag is set
       alns[0].aux(OFF_TARGET).is_ok() {    // skip if off-target flag is set (redundant with READ_FAILURE)
        return Ok(());
    }

    // skip reads in untargeted samples that map to other than nuclear chromosomes
    if let Some(writer) = writers.get_mut(&(tid as u32)){

        // commit on-target reads to temporary BAM files
        ctrs.increment(N_USABLE_READS);
        if alns.len() > 1 { ctrs.increment(N_SV_READS); }
        ctrs.increment_keyed(N_READS_BY_SAMPLE, sample_name);
        for aln in alns.iter_mut() {

            // add the sample bit tag for multi-sample comparison
            aln.push_aux(SB_TAG, Aux::U32(sample_bit)).unwrap();

            // commit on-target reads to temporary BAM files
            writer.write(aln)?; 
        }
        *read_counts.get_mut(&(tid as u32)).unwrap() += 1; // count reads, not alns

        // collect deduplicated alignment counts
        // rather than a CPU-intensive HashMap, uses a memory-intensive Vec sorted and chunked below
        let (unique_insert_span, was_reordered) = 
            UniqueInsertSpan::from_bam_record(&alns[0], sample_bit);
        let n_alns = alns.len();
        for i in 0..n_alns {
            let cigar_view = alns[i].cigar();
            channel_alns.push((
                ChannelAlignment {
                    unique_insert_span,
                    aln_i: (if was_reordered { n_alns - 1 - i } else { i }) as u8,
                },
                cigar_view.end_pos() as u32 - alns[i].pos() as u32
            ));
        }
    }
    Ok(())
}

// calculate on-target read and base counts over a (deduplicated) library
fn tally_read_bases(
    channel_alns: &mut Vec<(ChannelAlignment, u32)>,
    chroms:       &Chroms,
    w:            &mut Workflow,
    samples:      FxHashMap<u32, String>,
){
    let deduplicate_reads   = *w.cfg.get_bool(DEDUPLICATE_READS);
    let is_composite_genome = *w.cfg.get_bool(IS_COMPOSITE_GENOME);

    w.log.print("sorting alignments for deduplication");
    channel_alns.par_sort_unstable_by(|a, b| a.0.cmp(&b.0));

    w.log.print("tallying coverage statistics");
    for chunk in channel_alns.chunk_by(|a, b| a.0 == b.0) {
        let sample_name = samples.get(&chunk[0].0.unique_insert_span.sample_bit).unwrap();

        // get deduplicated counts
        let (n_dedup_instances, n_dedup_bases) = get_chunk_vals(chunk, deduplicate_reads);
        w.ctrs.increment(N_UNIQ_ALNS);
        w.ctrs.add_to(N_ALNS,  n_dedup_instances);
        w.ctrs.add_to(N_BASES, n_dedup_bases);
        w.ctrs.add_to_keyed(N_BASES_BY_SAMPLE,  sample_name, n_dedup_bases);

        // count reads and bases per genome to determine the observed mixing ratio
        if is_composite_genome {
            let channel_aln = &chunk[0].0;
            let ( chrom1, chrom_index1, _, _) = 
                SamRecord::unpack_signed_node(channel_aln.unique_insert_span.node1, &chroms);
            let (_chrom2, chrom_index2, _, _) = 
                SamRecord::unpack_signed_node(channel_aln.unique_insert_span.node2, &chroms);
            if chrom_index1 == chrom_index2 { // only count reads whose outer nodes map to the same genome
                let genome_name = chrom1.split_once('_').unwrap().1; // e.g., chr1_hs1
                if channel_aln.aln_i == 0 { // count reads, not alns
                    w.ctrs.increment_keyed(N_READS_BY_GENOME, genome_name);
                }
                w.ctrs.add_to_keyed(N_BASES_BY_GENOME,  genome_name, n_dedup_bases);
            }
        }
    }
}

// get (deduplicated) counts for a chunk of identical ChannelAlignments
fn get_chunk_vals(
    chunk: &[(ChannelAlignment, u32)],
    deduplicate_reads: bool,
) -> (usize, usize){ // n_dedup_instances, n_dedup_bases
    let n_instances = chunk.len();
    if deduplicate_reads || n_instances == 1 {
        (1, chunk[0].1 as usize)
    } else {
        let mut n_bases = 0_u32;
        for (_channel_aln, n_bases_aln) in chunk {
            n_bases += *n_bases_aln;
        }
        (n_instances, n_bases as usize)
    }
}

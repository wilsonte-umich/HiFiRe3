//! Split an input name-sorted BAM file into temporary per-chromosome BAM 
//! files based on the chromosome of the first alignment in each read.
//! 
//! Output only includes usable on-target reads, as defined by the absence 
//! the UNMAPPED flag bit and the READ_IS_OFF_TARGET and READ_FAILURE_FLAG tags.
//! 
//! Along the way, collect de-duplicated on-target coverage statistics.

// dependencies
use std::error::Error;
use std::fs::File;
use std::io::Write;
use rustc_hash::FxHashMap;
use rust_htslib::bam::{Reader, Read, Writer, Record as BamRecord, Header, Format};
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use genomex::genome::{Chroms, TargetRegions};
use genomex::sam::SamRecord;
use crate::formats::hf3_tags::*;
use crate::inserts::{UniqueInsertSpan, ChannelAlignment};

// constants
const TOOL: &str = "split_by_chrom";
pub_key_constants!(
    // from environment variables
    DEDUPLICATE_READS
    IS_COMPOSITE_GENOME
    NAME_BAM_FILE
    INDEX_FILE_PREFIX_WRK
    // counter keys
    N_READS
    N_USABLE_READS
    N_UNIQ_ALNS
    N_ALNS
    N_REF_BASES
    N_READ_BASES
    N_READS_BY_GENOME
    N_REF_BASES_BY_GENOME
    N_READ_BASES_BY_GENOME
);
const READ_FAILURE: &[u8] = READ_FAILURE_FLAG.as_bytes();
const OFF_TARGET: &[u8]   = READ_IS_OFF_TARGET.as_bytes(); 

/// ChannelAlignmentCounts holds (deduplicated) tallies for a ChannelAlignment. The 
/// value is a count or sum when DEDUPLICATE_READS is false, otherwise n_alns==1  
/// and n_ref_bases and n_read_bases are for the first encountered alignment.
struct ChannelAlignmentCounts {
    n_alns:      usize,
    n_ref_bases:  usize,
    n_read_bases: usize,
}

// main function called by xxx_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_bool_env(&[DEDUPLICATE_READS, IS_COMPOSITE_GENOME]);
    cfg.set_string_env(&[NAME_BAM_FILE, INDEX_FILE_PREFIX_WRK]);

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_READS,           "reads processed"),
        (N_USABLE_READS, "usable on-target reads in output"),
        (N_UNIQ_ALNS,       "unique (deduplicated) alignments in on-target reads"),
        (N_ALNS,            "total  (deduplicated) alignments in on-target reads"),
        (N_REF_BASES,       "(deduplicated) reference bases in on-target alignments"),
        (N_READ_BASES,      "(deduplicated) read bases in on-target alignments"),
    ]);
    ctrs.add_keyed_counters(&[
        (N_READS_BY_GENOME,      "on-target reads by genome"),
        (N_REF_BASES_BY_GENOME,  "reference bases in on-target alignments by genome"),
        (N_READ_BASES_BY_GENOME, "read bases in on-target alignments by genome"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();

    // collect the working chromosomes
    let chroms = Chroms::new(&mut w.cfg);
    let targets = TargetRegions::from_env(&mut w, false);
    let on_target_chroms = targets.get_region_chroms(&chroms);

    // initialize the output BAM writers, one per target chromosome
    let name_bam_path     = w.cfg.get_string(NAME_BAM_FILE);
    let chrom_file_prefix = w.cfg.get_string(INDEX_FILE_PREFIX_WRK);
    let mut name_bam = Reader::from_path(&name_bam_path)?;
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
        writers.insert(tid, Writer::from_path(
            format!("{}.chr{}.bam", chrom_file_prefix, chrom_index_padded),
            &header,
            Format::Bam
        ).expect(&format!("Failed to create BAM writer for chrom {}", chrom)));
        read_counts.insert(tid, 0); 
        chrom_map.insert(tid, chrom_index_padded);
    }

    // initialize the deduplication counter
    let deduplicating = *w.cfg.get_bool(DEDUPLICATE_READS);
    let mut insert_counts: FxHashMap<ChannelAlignment, ChannelAlignmentCounts> = FxHashMap::default();

    // process input BAM records
    w.log.print("streaming BAM records");
    let mut records = name_bam.records();
    let mut alns: Vec<BamRecord> = vec![records.next().unwrap()?]; // expect there to always be data, thus alns is never empty
    for result in records {
        let record = result?;
        if record.qname() != alns[0].qname() { 
            print_alns(&alns, &mut writers, &mut read_counts, &mut w.ctrs, deduplicating, &mut insert_counts)?;
            alns.clear();
        }
        alns.push(record);
    }
    print_alns(&alns, &mut writers, &mut read_counts, &mut w.ctrs, deduplicating, &mut insert_counts)?;

    // write chrom read count files, for setting capacity in downstream tools
    w.log.print("writing read count files");
    for (tid, count) in &read_counts {
        let chrom_index_padded = chrom_map.get(tid).unwrap();
        let count_file_path = format!("{}.chr{}.read_count", chrom_file_prefix, chrom_index_padded);
        let mut count_file = File::create(&count_file_path)?;
        writeln!(count_file, "{}", count)?;
    }
    
    // report counter values
    w.log.print("tallying coverage statistics");
    tally_read_bases(&insert_counts, &chroms, &mut w);
    w.ctrs.print_grouped(&[
        &[N_READS, N_USABLE_READS],
        &[N_UNIQ_ALNS, N_ALNS, N_REF_BASES, N_READ_BASES],
        &[N_READS_BY_GENOME],
        &[N_REF_BASES_BY_GENOME],
        &[N_READ_BASES_BY_GENOME],
    ]);
    Ok(())
}

// print and count on-target alignments
fn print_alns(
    alns:          &[BamRecord], 
    writers:       &mut FxHashMap<u32, Writer>,
    read_counts:   &mut FxHashMap<u32, usize>,
    ctrs:          &mut Counters,
    deduplicating: bool,
    insert_counts: &mut FxHashMap<ChannelAlignment, ChannelAlignmentCounts>
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
        for aln in alns { writer.write(aln)?; }
        *read_counts.get_mut(&(tid as u32)).unwrap() += 1; // count reads, not alns

        // collect deduplicated alignment counts
        let (unique_insert_span, was_reordered) = 
            UniqueInsertSpan::from_bam_record(&alns[0]);
        let n_alns = alns.len();
        for i in 0..n_alns {
            let channel_aln = ChannelAlignment {
                aln_i: (if was_reordered { n_alns - 1 - i } else { i }) as u8,
                unique_insert_span,
            };
            let counts = insert_counts.entry(channel_aln)
                .or_insert(ChannelAlignmentCounts {
                    n_alns:       0,
                    n_ref_bases:  0,
                    n_read_bases: 0,
                });
            if !deduplicating || counts.n_alns == 0 {
                let cigar_view = alns[i].cigar();
                counts.n_alns       += 1;
                counts.n_ref_bases  += cigar_view.end_pos() as usize - alns[i].pos() as usize;
                counts.n_read_bases += cigar_view.iter().fold(0usize, |acc, cig| {
                    match cig.char() {
                        'M' | '=' | 'X' | 'I' => acc + cig.len() as usize,
                        _ => acc,
                    }
                });
            }
        }
    }
    Ok(())
}

// calculate on-target read and base counts over a (deduplicated) library
fn tally_read_bases(
    insert_counts: &FxHashMap<ChannelAlignment, ChannelAlignmentCounts>,
    chroms: &Chroms,
    w: &mut Workflow,
){
    let is_composite_genome = *w.cfg.get_bool(IS_COMPOSITE_GENOME);
    for (channel_aln, counts) in insert_counts.iter() {
        w.ctrs.increment(N_UNIQ_ALNS);
        w.ctrs.add_to(N_ALNS,       counts.n_alns);
        w.ctrs.add_to(N_REF_BASES,  counts.n_ref_bases);
        w.ctrs.add_to(N_READ_BASES, counts.n_read_bases);

        // count reads and bases per genome to determine the observed mixing ratio
        if is_composite_genome {
            let ( chrom1, chrom_index1, _, _) = 
                SamRecord::unpack_signed_node(channel_aln.unique_insert_span.node1, &chroms);
            let (_chrom2, chrom_index2, _, _) = 
                SamRecord::unpack_signed_node(channel_aln.unique_insert_span.node2, &chroms);
            if chrom_index1 == chrom_index2 { // only count reads whose outer nodes map to the same genome
                let genome_name = chrom1.split_once('_').unwrap().1; // e.g., chr1_hs1
                if channel_aln.aln_i == 0 {
                    w.ctrs.increment_keyed(N_READS_BY_GENOME, genome_name);
                }
                w.ctrs.add_to_keyed(N_REF_BASES_BY_GENOME,  genome_name, counts.n_ref_bases);
                w.ctrs.add_to_keyed(N_READ_BASES_BY_GENOME, genome_name, counts.n_read_bases);
            }
        }
    }
}

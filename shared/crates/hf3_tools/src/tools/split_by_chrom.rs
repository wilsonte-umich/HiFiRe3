//! Split an input name-sorted BAM file into temporary per-chromosome BAM 
//! files based on the chromosome of the first alignment in each read.
//! 
//! Output only includes on-target reads, as defined by the absence of the
//! READ_IS_OFF_TARGET tag.
//! 
//! Along the way, collect de-duplicated on-target coverage statistics.

// dependencies
use std::error::Error;
use std::collections::HashMap;
use rust_htslib::bam::{Reader, Read, Writer, Record, Header, Format, record::Aux};
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use genomex::genome::{Chroms, TargetRegions};
use genomex::sam::SamRecord;
use genomex::bam::tags;
use crate::formats::hf3_tags::*;

// constants
const TOOL: &str = "split_by_chrom";
pub_key_constants!(
    // from environment variables
    NAME_BAM_FILE
    INDEX_FILE_PREFIX_WRK
    DEDUPLICATE_READS
    IS_COMPOSITE_GENOME
    // counter keys
    N_READS
    N_READS_ON_TARGET
    N_UNIQ_ALNS
    N_ALNS
    N_REF_BASES
    N_READ_BASES
    N_READS_BY_GENOME
    N_REF_BASES_BY_GENOME
    N_READ_BASES_BY_GENOME
);
const OFF_TARGET_FLAG: &[u8] = READ_IS_OFF_TARGET.as_bytes(); 

/// ChannelAlignment holds (re)ordered alignment indices, (re)ordered outer nodes,
/// and ONT channel information for deduplicated counting.
/// 
/// Thus, an instance identifies "the nth alignment of a read after (re)ordering into the 
/// canonical orientation of the outer nodes (for an insert on a specific ONT channel)".
#[derive(Hash, Eq, PartialEq)]
struct ChannelAlignment {
    alni:    u8, // up to 256 alignments per read
    node1:   isize,
    node2:   isize,
    channel: u32,
}
/// ChannelAlignmentCounts holds (deduplicated) tallies for a ChannelAlignment. The 
/// value is a count or sum when DEDUPLICATE_READS is false, otherwise n==1 and 
/// n_ref_bases and n_read_bases are for the first encountered alignment.
struct ChannelAlignmentCounts {
    n:            usize,
    n_ref_bases:  usize,
    n_read_bases: usize,
}

// main function called by xxx_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[NAME_BAM_FILE, INDEX_FILE_PREFIX_WRK]);
    cfg.set_bool_env(&[DEDUPLICATE_READS, IS_COMPOSITE_GENOME]);

    // initialize counters
    let mut ctrs = Counters::new(TOOL, &[
        (N_READS,           "reads processed"),
        (N_READS_ON_TARGET, "on-target reads in output"),
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
    w.log.print("initializing");

    // collect the working chromosomes
    let chroms = Chroms::new(&mut w.cfg);
    let targets = TargetRegions::new(&mut w);
    let on_target_chroms = targets.get_on_target_chroms(&chroms);

    // initialize the output BAM writers, one per target chromosome
    let name_bam_path = w.cfg.get_string(NAME_BAM_FILE);
    let chrom_file_prefix = w.cfg.get_string(INDEX_FILE_PREFIX_WRK);
    let mut name_bam = Reader::from_path(&name_bam_path)?;
    let header_view = name_bam.header(); // for TID lookups
    let header = Header::from_template(header_view); // shared header for each output BAM writer
    let mut writers: HashMap<u32, Writer> = HashMap::new(); // keyed by TID, file named by our padded chrom index
    for (chrom, chrom_index) in on_target_chroms {
        let tid = header_view.tid(chrom.as_bytes()).unwrap();
        let chrom_index_padded = format!("{:02}", chrom_index);
        writers.insert(tid, Writer::from_path(
            format!("{}.chr{}.bam", chrom_file_prefix, chrom_index_padded),
            &header,
            Format::Bam
        ).expect(&format!("Failed to create BAM writer for chrom {}", chrom)));
    }

    // initialize the deduplication counter
    let deduplicating = *w.cfg.get_bool(DEDUPLICATE_READS);
    let mut insert_counts: HashMap<ChannelAlignment, ChannelAlignmentCounts> = HashMap::new();

    // process input BAM records
    let mut records = name_bam.records();
    w.log.print("streaming BAM records");
    let mut alns: Vec<Record> = vec![records.next().unwrap()?]; // expect there to always be data, thus alns is never empty for print_alns
    for result in records {
        let record = result?;
        if record.qname() != alns[0].qname() { 
            print_alns(&alns, &mut writers, &mut w.ctrs, deduplicating, &mut insert_counts)?;
            alns.clear();
        }
        alns.push(record);
    }
    print_alns(&alns, &mut writers, &mut w.ctrs, deduplicating, &mut insert_counts)?;

    // report counter values
    w.log.print("tallying coverage statistics");
    tally_read_bases(&insert_counts, &chroms, &mut w);
    w.ctrs.print_grouped(&[
        &[N_READS, N_READS_ON_TARGET],
        &[N_UNIQ_ALNS, N_ALNS, N_REF_BASES, N_READ_BASES],
        &[N_READS_BY_GENOME],
        &[N_REF_BASES_BY_GENOME],
        &[N_READ_BASES_BY_GENOME],
    ]);
    Ok(())
}

// print and count on-target alignments
fn print_alns(
    alns:          &[Record], 
    writers:       &mut HashMap<u32, Writer>,
    ctrs:          &mut Counters,
    deduplicating: bool,
    insert_counts: &mut HashMap<ChannelAlignment, ChannelAlignmentCounts>
) -> Result<(), Box<dyn Error>> {

    // skip off-target reads
    ctrs.increment(N_READS);
    let tid = alns[0].tid();
    if tid < 0 ||                                  // skip unmapped reads
       alns[0].aux(OFF_TARGET_FLAG).is_ok() { // skip if off-target flag is set
        return Ok(());
    }

    // commit on-target reads to temporary BAM files
    ctrs.increment(N_READS_ON_TARGET);
    let writer = writers.get_mut(&(tid as u32)).unwrap();
    for aln in alns {
        writer.write(aln)?;
    }

    // collect deduplication statistics
    let channel = tags::get_tag_u32_default(&alns[0], CHANNEL, 0);
    let ordered_outer_nodes = SamRecord::sam_tag_to_paired_nodes(
        &tags::get_tag_str(&alns[0], OUTER_NODES), 
        true
    );
    let n_alns = alns.len();
    for i in 0..n_alns {
        let channel_aln = ChannelAlignment {
            alni: (if ordered_outer_nodes.was_reordered { n_alns - 1 - i } else { i }) as u8,
            node1: ordered_outer_nodes.node1, 
            node2: ordered_outer_nodes.node2, 
            channel,
        };
        let counts = insert_counts.entry(channel_aln)
            .or_insert(ChannelAlignmentCounts {
                n:            0,
                n_ref_bases:  0,
                n_read_bases: 0,
            });
        if !deduplicating || counts.n == 0 {
            let cigar_view = alns[i].cigar();
            counts.n += 1;
            counts.n_ref_bases  += cigar_view.end_pos() as usize - alns[i].pos() as usize;
            counts.n_read_bases += alns[i].seq_len() - (cigar_view.leading_softclips() + cigar_view.trailing_softclips()) as usize;
        }
    }
    Ok(())
}

// calculate on-target read and base counts over a (deduplicated) library
fn tally_read_bases(
    insert_counts: &HashMap<ChannelAlignment, ChannelAlignmentCounts>,
    chroms: &Chroms,
    w: &mut Workflow,
){
    let is_composite_genome = *w.cfg.get_bool(IS_COMPOSITE_GENOME);
    for (channel_aln, counts) in insert_counts.iter() {
        w.ctrs.increment(N_UNIQ_ALNS);
        w.ctrs.add_to(N_ALNS, counts.n);
        w.ctrs.add_to(N_REF_BASES, counts.n_ref_bases);
        w.ctrs.add_to(N_READ_BASES, counts.n_read_bases);

        // count reads and bases per genome to determine the observed mixing ratio
        if is_composite_genome {
            let (chrom1, _, _) = SamRecord::unpack_signed_node(channel_aln.node1, &chroms);
            let (chrom2, _, _) = SamRecord::unpack_signed_node(channel_aln.node2, &chroms);
            if chrom1 == chrom2 {
                let genome_name = chrom1.split_once('_').unwrap().1; // e.g., chr1_hs1
                if channel_aln.alni == 0 {
                    w.ctrs.increment_keyed(N_READS_BY_GENOME, genome_name);
                }
                w.ctrs.add_to_keyed(N_REF_BASES_BY_GENOME,  genome_name, counts.n_ref_bases);
                w.ctrs.add_to_keyed(N_READ_BASES_BY_GENOME, genome_name, counts.n_read_bases);
            }
        }
    }
}

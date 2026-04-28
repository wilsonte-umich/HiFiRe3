//! Perform the pre- and post-alignment processing steps on Ultima aligned reads
//! that are performed on non-Ultima reads by prepare_fastq.pl, fastp, and minimap2.
//!    - tabulate observed insert size distributions, i.e., non-fixed read lengths
//!    - enforce PLATFORM_MIN_INSERT_SIZE
//!    - perform read quality filtering
//! All tasks can be performed in an alignment autonomous manner, although are
//! done more than once on reads with split alignments.

// dependencies
use std::error::Error;
use rustc_hash::FxHashMap;
use mdi::OutputFile;
use rust_htslib::bam::{Reader, Read, Record, Writer, Header, Format, CompressionLevel};
use rust_htslib::tpool::ThreadPool;
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};

// constants
const TOOL: &str = "post_process_ultima";
pub_key_constants!(
    // from environment variables
    N_CPU
    READ_1_FILES
    GENOME_FASTA
    QUALIFIED_QUALITY_PHRED
    UNQUALIFIED_PERCENT_LIMIT
    AVERAGE_QUAL
    PLATFORM_MIN_INSERT_SIZE
    UNFILTERED_INSERT_SIZES_FILE
    // counter keys
    N_ALNS_IN
    N_ALNS_OUT
);
pub const SIZE_PLOT_BIN_SIZE: usize = 10;

/// Tool carrying processing options.
struct Tool {
    qualified_quality_phred:    u8,
    unqualified_fraction_limit: f64,
    average_qual:               f64,
    platform_min_insert_size:   usize,
    bam_out:                    Writer,
    insert_size_counts:         FxHashMap<usize, usize>,
}

// main function called by xxx_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_u8_env(&[QUALIFIED_QUALITY_PHRED, UNQUALIFIED_PERCENT_LIMIT, AVERAGE_QUAL]);
    cfg.set_u32_env(&[N_CPU]);
    cfg.set_usize_env(&[PLATFORM_MIN_INSERT_SIZE]);
    cfg.set_string_env(&[READ_1_FILES, GENOME_FASTA]);

    // initialize counters
    let ctrs = Counters::new(TOOL, &[
        (N_ALNS_IN,  "alignments processed"),
        (N_ALNS_OUT, "alignments in output"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();

    // use a thread pool for CRAM reading and bam writing
    let tpool = ThreadPool::new(w.cfg.get_u32(N_CPU) - 1)?;

    // open the shared output writing to concatenate all input files
    let coord_cram_files = w.cfg.get_string(READ_1_FILES);
    let coord_cram_paths = coord_cram_files.split(',').collect::<Vec<&str>>();
    let coord_cram = Reader::from_path(coord_cram_paths[0])?;
    let header = Header::from_template(coord_cram.header());
    let mut bam_out = Writer::from_stdout(&header, Format::Bam)?;
    bam_out.set_compression_level(CompressionLevel::Uncompressed)?;
    bam_out.set_thread_pool(&tpool)?;

    // get ready for alignment processing
    drop(coord_cram);
    let mut tool = Tool {
        qualified_quality_phred:    *w.cfg.get_u8(QUALIFIED_QUALITY_PHRED),
        unqualified_fraction_limit: (*w.cfg.get_u8(UNQUALIFIED_PERCENT_LIMIT) as f64) / 100.0,
        average_qual:               *w.cfg.get_u8(AVERAGE_QUAL) as f64,
        platform_min_insert_size:   *w.cfg.get_usize(PLATFORM_MIN_INSERT_SIZE),
        bam_out:                    bam_out,
        insert_size_counts:         FxHashMap::default(),
    };

    // sequentially process aligned input CRAM files
    let genome_fasta = w.cfg.get_string(GENOME_FASTA);
    for coord_cram_path in coord_cram_paths {
        let cram_file_name = coord_cram_path.split('/').last().unwrap();
        eprintln!("    {}", cram_file_name);
        let mut coord_cram = Reader::from_path(coord_cram_path)?;
        coord_cram.set_reference(genome_fasta)?;
        coord_cram.set_thread_pool(&tpool)?;
        
        // process alignments one at a time
        let mut aln = Record::new();
        while let Some(result) = coord_cram.read(&mut aln) {
            match result {
                Ok(_) => {
                    process_aln(&aln, &mut tool, &mut w.ctrs)?;
                }
                Err(_) => panic!("BAM parsing failed.")
            }
        }
    }

    // print unfiltered insert size distribution to file for plotting
    let mut insert_sizes_file = OutputFile::open_env(&mut w.cfg, UNFILTERED_INSERT_SIZES_FILE);

    let mut insert_sizes = tool.insert_size_counts.keys().cloned().collect::<Vec<usize>>();
    insert_sizes.sort_unstable();
    for insert_size in insert_sizes {
        let count = tool.insert_size_counts.get(&insert_size).unwrap();
        insert_sizes_file.write_record(vec![
            "unfiltered",
            &insert_size.to_string(),
            &count.to_string(),
        ]);
    }

    // report counter values
    w.ctrs.print_grouped(&[
        &[N_ALNS_IN, N_ALNS_OUT],
    ]);
    Ok(())
}

fn process_aln(
    aln:  &Record,
    tool: &mut Tool,
    ctrs: &mut Counters
) -> Result<(), Box<dyn Error>> {
    ctrs.increment(N_ALNS_IN);

    // record all read insert sizes after upstream adapter trimming but before quality filtering
    let seq_len = aln.seq_len();
    if !aln.is_supplementary() {
        let size_bin = (seq_len / SIZE_PLOT_BIN_SIZE) * SIZE_PLOT_BIN_SIZE;
        *tool.insert_size_counts.entry(size_bin).or_insert(0) += 1;
    }

    // enforce min insert size
    if seq_len < tool.platform_min_insert_size { 
        return Ok(());
    }
    let seq_len_f64 = seq_len as f64;

    // enforce min quality metrics
    let mut qual_sum: f64 = 0.0;
    let mut n_bases_failed_qual: f64 = 0.0;
    aln.qual().iter().for_each(|q|{
        qual_sum += *q as f64;
        if *q < tool.qualified_quality_phred {
            n_bases_failed_qual += 1.0;
        }
    });
    if qual_sum / seq_len_f64 < tool.average_qual ||
       n_bases_failed_qual > seq_len_f64 * tool.unqualified_fraction_limit {
        return Ok(());
    }

    // commit kept alns
    tool.bam_out.write(aln)?;
    ctrs.increment(N_ALNS_OUT);
    Ok(())
}

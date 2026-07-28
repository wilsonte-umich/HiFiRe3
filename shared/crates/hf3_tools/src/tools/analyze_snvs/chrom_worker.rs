//! Process reads with first alignments on a specific chromosome, provided as a 
//! message on a channel.

// imports
use std::error::Error;
use crossbeam::channel::{Receiver, Sender};
use minimap2::{Aligner as Minimap2};
use rust_htslib::bam::{Reader, Read, Record as BamRecord};
use mdi::pub_key_constants;
use mdi::workflow::Config;
use crate::sites::SiteMatches;
use crate::snvs::*;

// constants
pub_key_constants!(
    // from environment variables
    INDEX_FILE_PREFIX_WRK
    GENOME_REPEAT_MASKER_BED
    GENOME_SIMPLE_REPEAT_BED
);
const MAX_CLIP: u32 = 25;

// process chromosomes received on the channel
pub fn process_chrom(
    tool:     &SnvAnalysisTool,
    rx_chrom: Receiver<(String, u8)>,
    tx_data:  Sender<SnvChromWorkerData>,
) -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[INDEX_FILE_PREFIX_WRK, GENOME_REPEAT_MASKER_BED, GENOME_SIMPLE_REPEAT_BED]);
    let chrom_file_prefix = cfg.get_string(INDEX_FILE_PREFIX_WRK); // created by split_bam_by_chrom
    let rmsk_simple_repeats_bed = cfg.get_string(GENOME_REPEAT_MASKER_BED);
    let trf_simple_repeats_bed = cfg.get_string(GENOME_SIMPLE_REPEAT_BED);

    // process chromosomes received on the channel
    for (chrom_name, chrom_index) in rx_chrom.iter() {
        let chrom_index_padded = format!("{:02}", chrom_index);
        let chrom_bam_path  = format!("{}.chr{}.bam", chrom_file_prefix, &chrom_index_padded);
        eprintln!("    {}", chrom_name);

        // open the input BAM file
        // all reads are on-target and have first alignment on chrom
        let mut chrom_bam = Reader::from_path(&chrom_bam_path)?;

        // assemble the chromosome worker tool
        let variants_file_path = format!(
            "{}.chr{}.snv_indel.txt.bgz", 
            chrom_file_prefix, &chrom_index_padded
        );
        let variant_reads_file_path = format!(
            "{}.chr{}.variant_reads.txt.bgz", 
            chrom_file_prefix, &chrom_index_padded
        );
        let reads_on_reference_path = format!(
            "{}.chr{}.encodings.reads_on_reference.bed.bgz", 
            chrom_file_prefix, &chrom_index_padded
        );
        let reads_on_haplotype_path = format!(
            "{}.chr{}.encodings.reads_on_haplotype.bed.bgz", 
            chrom_file_prefix, &chrom_index_padded
        );
        let mut worker = SnvChromWorker {
            chrom: chrom_name.clone(),
            chrom_index,
            chrom_tid: tool.fa.fai().tid(&chrom_name).expect("Failed to get chrom TID"),
            simple_repeats: SimpleRepeats::new(
                tool, &chrom_name, rmsk_simple_repeats_bed, trf_simple_repeats_bed
            ),
            minimap2: Minimap2::builder().map_hifi().with_cigar(),
            variant_tally:       VariantsTally::new(),
            variant_reads_tally: VariantReadsTally::new(),
            reads_on_reference:  AlignmentEncodings::new(),
            reads_on_haplotype:  AlignmentEncodings::new(),
            variants_file_path,
            variant_reads_file_path,
            reads_on_reference_path,
            reads_on_haplotype_path,
        };

        // process alignment records one at a time, add to growing RE fragment collections
        let mut aln = BamRecord::new();
        let mut chrom_aln_count:    usize = 0;
        let mut chrom_aln_count_used: usize = 0;
        let mut fragment_reads = FragmentReads::new();
        while let Some(result) = chrom_bam.read(&mut aln) {
            match result {
                Ok(_)  => {
                    chrom_aln_count += 1;
                    process_aln(&aln, &mut chrom_aln_count_used, &mut fragment_reads)?;
                },
                Err(_) => panic!("BAM parsing failed")
            }
        }

        // post-process read groups by re-aligning reads to fragment consensus(es)
        let mut haplotype_consensuses = HaplotypeConsensuses::new();
        worker.analyze_reads(tool, fragment_reads, &mut haplotype_consensuses);

        // finish processing and writing pileup and variants
        let variant_metadata = VariantsTally::write_sorted(tool, &mut worker);
        let variant_reads_metadata = VariantReadsTally::write_sorted(tool, &mut worker);
        let clonal_metadata = AlignmentEncodings::write_sorted(
            tool, &worker, &mut haplotype_consensuses,
            &worker.reads_on_reference, 
            &worker.reads_on_reference_path
        );
        let subclonal_metadata = AlignmentEncodings::write_sorted(
            tool, &worker, &mut haplotype_consensuses,
            &worker.reads_on_haplotype, 
            &worker.reads_on_haplotype_path
        );

        // send error corrected metadata to main thread
        tx_data.send(SnvChromWorkerData::TotalAlnCount(chrom_aln_count))?;
        tx_data.send(SnvChromWorkerData::UsableAlnCount((chrom_name.clone(), chrom_aln_count_used)))?;
        tx_data.send(SnvChromWorkerData::VariantMetadata(variant_metadata))?;
        tx_data.send(SnvChromWorkerData::VariantReadsMetadata(variant_reads_metadata))?;
        tx_data.send(SnvChromWorkerData::ClonalEncodingMetadata(clonal_metadata))?;
        tx_data.send(SnvChromWorkerData::SubclonalEncodingMetadata(subclonal_metadata))?;
    }
    Ok(())
}

// process one alignment
fn process_aln(
    aln: &BamRecord, 
    chrom_aln_count_used: &mut usize,
    fragment_reads: &mut FragmentReads,
) -> Result<(), Box<dyn Error>>{

    // short-circuit rare reads that do not yield a productive RE fragment match
    let sites = SiteMatches::from_bam_record(&aln);
    if sites.site5.pos1 == 0 || 
       sites.site3.pos1 == 0 || 
       sites.site5.pos1 == sites.site3.pos1 || 
       sites.site5.distance.unsigned_abs() > MAX_CLIP || 
       sites.site3.distance.unsigned_abs() > MAX_CLIP {
        return Ok(());
    }
    *chrom_aln_count_used += 1;

    // collect reads by RE fragment for later dismissal or alignment processing
    let re_fragment = ReFragment::new(
        sites.site5.pos1, sites.site3.pos1, aln.is_reverse()
    );
    let read_instance = ReadInstance::new(aln);
    fragment_reads.insert(re_fragment, read_instance);

    Ok(())
}

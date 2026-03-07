//! Process reads with first alignments on a specific chromosome, 
//! provided as a message on a channel.

// dependencies
use std::error::Error;
use crossbeam::channel::{Receiver, Sender};
use rust_htslib::bam::{Reader, Read, Record as BamRecord};
use mdi::pub_key_constants;
use mdi::workflow::Config;
use genomex::bam::tags as bam_tags;
use genomex::sam::cigar::CigarString;
use crate::formats::hf3_tags::*;
use crate::snvs::{
    tags as snv_tags, 
    SnvChromWorkerData, SnvChromWorker, 
    SnvAnalysisTool, 
    ChromPileup, ChromVariants
};

// constants
pub_key_constants!(
    // from environment variables
    INDEX_FILE_PREFIX_WRK
    MIN_N_PASSES
    MIN_SNV_INDEL_QUAL
);

// process chromosomes received on the channel
pub fn process_chrom(
    tool:     &SnvAnalysisTool,
    rx_chrom: Receiver<(String, u8)>,
    tx_data:  Sender<SnvChromWorkerData>,
) -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_u8_env(&[MIN_N_PASSES, MIN_SNV_INDEL_QUAL]);
    cfg.set_string_env(&[INDEX_FILE_PREFIX_WRK]);
    let chrom_file_prefix = cfg.get_string(INDEX_FILE_PREFIX_WRK);

    // process chromosomes received on the channel
    for (chrom_name, chrom_index) in rx_chrom.iter() {
        let chrom_index_padded = format!("{:02}", chrom_index);
        let chrom_size = *tool.chroms.sizes.get(&chrom_name).unwrap();

        // run SNV/indel analysis twice
        //   once with all reads for maximally sensitive (sub)clonal variant calling
        //   once with only error-corrected reads filtered for n_passes for maximally specific rare variant calling
        for read_type in ["all_reads", "error_corrected"] {
            eprintln!("    {} {}", chrom_name, read_type);

            // open the input BAM file
            // all reads are on-target and have first alignment on chrom
            let chrom_bam_path  = format!("{}.chr{}.bam", chrom_file_prefix, &chrom_index_padded);
            let mut chrom_bam = Reader::from_path(&chrom_bam_path)?;

            // assemble the chromosome worker tool
            let pileup_file_path   = format!(
                "{}.chr{}.{}.pileup.bed.bgz",    
                chrom_file_prefix, &chrom_index_padded, read_type
            );
            let variants_file_path = format!(
                "{}.chr{}.{}.snv_indel.txt.bgz", 
                chrom_file_prefix, &chrom_index_padded, read_type
            );
            let mut worker = SnvChromWorker {
                chrom:              chrom_name.clone(),
                chrom_index,
                include_all_reads:  read_type == "all_reads",
                min_n_passes:       *cfg.get_u8(MIN_N_PASSES),
                min_snv_indel_qual: *cfg.get_u8(MIN_SNV_INDEL_QUAL) as usize,
                pileup:             ChromPileup::with_capacity(chrom_size as usize),
                variants:           ChromVariants::new(),
                pileup_file_path,
                variants_file_path,
            };

            // process alignment records one at a time
            let mut aln = BamRecord::new();
            let mut chrom_aln_count:    usize = 0;
            let mut chrom_aln_count_ec: usize = 0;
            while let Some(result) = chrom_bam.read(&mut aln) {
                match result {
                    Ok(_)  => {
                        process_aln(&aln, &mut worker, &mut chrom_aln_count_ec)?;
                        chrom_aln_count += 1;
                    },
                    Err(_) => panic!("BAM parsing failed")
                }
            }

            // finish processing and writing pileup and variants
            let pileup_metadata   = ChromPileup::write_chunked( tool, &mut worker);
            let variant_metadata = ChromVariants::write_sorted(tool, &mut worker);

            // send chrom read count to main thread
            if !worker.include_all_reads {
                tx_data.send(SnvChromWorkerData::TotalAlnCount(chrom_aln_count))?;
                tx_data.send(SnvChromWorkerData::ErrorCorrectedAlignmentCount((chrom_name.clone(), chrom_aln_count_ec)))?;
                tx_data.send(SnvChromWorkerData::PileupMetadata(pileup_metadata))?;
                tx_data.send(SnvChromWorkerData::VariantMetadata(variant_metadata))?;
            }
        }
    }
    Ok(())
}

// process one alignment
fn process_aln(
    aln:    &BamRecord, 
    worker: &mut SnvChromWorker,
    chrom_aln_count_ec: &mut usize,
) -> Result<(), Box<dyn Error>>{

    // short-circuit reads that do not meet duplex filtering criteria when restricting to error-corrected reads
    // valid EC reads must have a dd tag and meet the minimum number of total PacBio passes reqested by the user
    let dd = bam_tags::get_tag_str_opt(aln, STRAND_DIFFERENCES);
    let n_passes = bam_tags::get_tag_f32_default(aln, PACBIO_EFF_COVERAGE, 0.0) as u8;
    if !worker.include_all_reads {
        if dd.is_none() || 
           n_passes < worker.min_n_passes {
            return Ok(());
        } else {
            *chrom_aln_count_ec += 1;
        }
    }

    // when available, process the read-level dd:Z: tag into a variant calling mask
    let cigar      = aln.cigar().to_string();
    let cigar = CigarString{cigar: cigar};
    let read_len    = cigar.get_read_len() as usize;
    let mask = if let Some(dd) = dd {
        let is_reverse = aln.is_reverse(); // expect most but not all alignments with a dd tag to be forward
        snv_tags::DdTag(dd).get_read_mask(read_len, is_reverse)   // in alignment order, not read order
    } else {
        vec![false; read_len]
    };

    // process the alignment-level cs:Z: tag into a read pileup and allowed variant list
    let left_clip = cigar.get_clip_left();
    let ref_pos0 = aln.pos() as u32;
    let cs = bam_tags::get_tag_str(aln, DIFFERENCE_STRING);
    let cs = snv_tags::CsTag(cs);
    let sample_bit = bam_tags::get_tag_u16(aln, SAMPLE_BIT);
    cs.process_aln(
        worker, 
        &mask, 
        aln.qual(),
        left_clip, 
        ref_pos0, 
        sample_bit, 
        n_passes
    ); 
    Ok(())
}

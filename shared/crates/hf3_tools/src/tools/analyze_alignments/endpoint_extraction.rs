//! Module to collect best endpoints of sequenced inserts in the 
//! reference genome as candidate RE site positions.

// dependencies
use std::error::Error;
use rustc_hash::FxHashMap;
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use mdi::OutputFile;
use genomex::sam::{SamRecord, flag};
use genomex::genome::Chroms;
use crate::formats::hf3_tags::TRIM_LENGTHS;

// constants
pub_key_constants!{
    // from environment variables
    EXPECTING_ENDPOINT_RE_SITES
    MIN_MAPQ
    CLIP_TOLERANCE
    // counter keys
    N_ENDS
    N_HIGH_QUAL
    N_ENDPOINTS
    // output files
    OBSERVED_ENDPOINTS_FILE
}
const ONT_ADAPTER_LEN5: i32 = 34; // 3 empirically determined median adapter trim lengths 
const TRIM5: usize = 0;           // tl tag indices for 5' and 3' adapter trims
const TRIM3: usize = 1;
const FORWARD: u8  = 0;
const REVERSE: u8  = 1;

/// Endpoint structure for collecting and printing putative RE site positions.
pub struct Endpoints {
    collecting:             bool,
    min_mapq:               u8,
    clip_tolerance:         u32,
    is_ont:                 bool,
    is_end_to_end_platform: bool,
    is_paired_reads:        bool,
    counter:                FxHashMap<(u8, u32), usize>, // (chrom_index, site_pos1) -> count
    file:                   OutputFile,
}
impl Endpoints {

    /// Initialize endpoint extraction, opening a file prior to processing 
    /// in case IO fails.
    pub fn new(w: &mut Workflow) -> Endpoints {
        w.cfg.set_u8_env(&[MIN_MAPQ]);
        w.cfg.set_u32_env(&[CLIP_TOLERANCE]);
        w.cfg.set_bool_env(&[EXPECTING_ENDPOINT_RE_SITES]);
        let collecting = *w.cfg.get_bool(EXPECTING_ENDPOINT_RE_SITES);
        if collecting {
            w.ctrs.add_counters(&[
                (N_ENDS, "read endpoints processed"),
                (N_HIGH_QUAL, format!(
                    "candidate endpoints with MAPQ >= {} and clip <= {}", 
                    w.cfg.get_u8(MIN_MAPQ), w.cfg.get_u32(CLIP_TOLERANCE)).as_str()
                ),
                (N_ENDPOINTS, "unique endpoint reference positions nominated as RE sites"),
            ]);
        }
        Endpoints {
            collecting,
            min_mapq: *w.cfg.get_u8(MIN_MAPQ),
            clip_tolerance: *w.cfg.get_u32(CLIP_TOLERANCE),
            is_ont: *w.cfg.get_bool(super::IS_ONT),
            is_end_to_end_platform: *w.cfg.get_bool(super::IS_END_TO_END_PLATFORM),
            is_paired_reads: *w.cfg.get_bool(super::IS_PAIRED_READS),
            counter: FxHashMap::default(), // (chrom_index, site_pos1) -> count, incompatible with mdi::keyed_counters
            file: OutputFile::open_env(&mut w.cfg, OBSERVED_ENDPOINTS_FILE),
        }
    }

    /// Process the alignments for one read (of a pair).
    pub fn extract(
        &mut self,
        read_n: usize, 
        alns: &mut [SamRecord], 
        w: &mut Workflow,
        chroms: &Chroms,
        is_unmerged_pair: bool,
    ) -> Result<(), Box<dyn Error>> {
        if !self.collecting { return Ok(()); }

        // get ONT adapter trims
        let trims: Vec<usize> = if self.is_ont {
            if let Some(tl) = alns[0].get_tag_value(TRIM_LENGTHS) {
                tl.split(',').map(|x| x.parse().unwrap()).collect()
            } else {
                vec![0, 0]
            }
        } else {
            vec![0, 0]
        };

        // always attempt to process the 5' end of all reads
        self.process_endpoint(&alns[0], &trims, 5, w, chroms)?;

        // the 3' end of paired reads is irrelevant
        if read_n == 2 { return Ok(()); }

        // attempt to process read 3' ends if informative, i.e., inferred to be from end-to-end read
        // NOT necessarily on the same chromosome or part of the same contig as the 5' end
        if self.is_end_to_end_platform || // PacBio or other platforms that guarantee end-to-end reads
            (
                self.is_ont && // complete ONT reads, must be trimmed at 3' end to be confident end-to-end
                trims[TRIM3] > 0
            ) ||
            (
                self.is_paired_reads &&  // merged paired reads
                !is_unmerged_pair        // the two ends of unmerged pairs both captured as 5' ends above
            ){
            self.process_endpoint(&alns[alns.len() - 1], &trims, 3, w, chroms)?;
        }
        Ok(())
    }

    // check an informative endpoint for sufficient alignment quality
    fn process_endpoint(
        &mut self,
        aln: &SamRecord, 
        trims: &Vec<usize>, 
        end: u8, 
        w: &mut Workflow,
        chroms: &Chroms,
    ) -> Result<(), Box<dyn Error>> {
        w.ctrs.increment(N_ENDS);

        // check for an informative endpoint with sufficient alignment quality
        // implicitly rejects unmapped with MAPQ=0
        if !chroms.is_canonical(&aln.rname) || aln.mapq < self.min_mapq { 
            return Ok(()); 
        }
        let strand = if aln.check_flag_any(flag::REVERSE) { REVERSE } else { FORWARD };
        let clip_len = if end == 5 { 
            aln.get_query_start0() 
        } else if strand == FORWARD { 
            aln.get_clip_right()
        } else { 
            aln.get_clip_left()
        };
        if clip_len > self.clip_tolerance { 
            return Ok(()); 
        }
        w.ctrs.increment(N_HIGH_QUAL);

        // parse the combination of endpoint rPos1, strand to predicted sitePos1
        //        *           *           *     * = sitePos1
        //   ----|xx5-----3xx/xx5-----3xx|----  as nodes are numbered for alignments
        //   ----|xx3-----5xx/xx3-----5xx|----  
        let mut site_pos1 = if end == 5 {
            if strand == FORWARD { aln.pos1 } else { aln.get_end1() + 1 }
        } else {
            if strand == FORWARD { aln.get_end1() + 1 } else { aln.pos1 }
        };

        // adjust sitePos1 based on terminal alignment clips
        let clip_dir = match (end, strand) {
            (5, FORWARD) => -1,
            (5, REVERSE) =>  1,
            (3, FORWARD) =>  1,
            (3, REVERSE) => -1,
            _ => unreachable!(),
        };
        site_pos1 = (site_pos1 as i32 + clip_dir * clip_len as i32) as u32; // concern for possible overflow

        // use a typical trim value to adjust clips on ONT trim failures, to account for adapter bases still present
        // only applicable to 5' ends here, untrimmed 3' ends were filtered above as uninformative
        if self.is_ont && end == 5 && trims[TRIM5] == 0 {
            site_pos1 = (site_pos1 as i32 - clip_dir * ONT_ADAPTER_LEN5) as u32;
        }

        // keep a tally of all genome positions nominated as RE sites
        *self.counter.entry((chroms.index[&aln.rname], site_pos1)).or_insert(0) += 1;
        Ok(())
    }

    /// Report all unique observed endpoint positions with counts.
    pub fn write(&mut self, w: &mut Workflow, chroms: &Chroms) -> Result<(), Box<dyn Error>> {
        if !self.collecting { return Ok(()); }
        w.log.print("printing endpoint tallies");
        let mut sorted_keys: Vec<(u8, u32)> = self.counter.keys().cloned().collect();
        sorted_keys.sort_unstable();
        for (chrom_index, site_pos1) in sorted_keys.iter() {
            w.ctrs.increment(N_ENDPOINTS);
            self.file.write_record(vec!(
                &chroms.rev_index[chrom_index],
                &site_pos1.to_string(),
                &self.counter[&(*chrom_index, *site_pos1)].to_string()
            ));
        }
        self.file.close();
        Ok(())
    }
}

//! Module to collect best endpoints of sequenced inserts in the 
//! reference genome as candidate RE site positions.

// dependencies
use std::error::Error;
use std::collections::HashMap;
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use mdi::OutputFile;
use genomex::sam::{SamRecord, flag};
use super::Tool;

// constants
pub_key_constants!{
    // from environment variables
    EXPECTING_ENDPOINT_RE_SITES
    MIN_MAPQ
    CLIP_TOLERANCE
    // derived values
    IS_END_TO_END_READ
    IS_PAIRED_READS
    IS_ONT
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
const FORWARD: u8 = 0;
const REVERSE: u8 = 1;

/// Endpoint structure for collecting and printing putative RE site positions.
pub struct Endpoints {
    collecting: bool,
    counter: HashMap<(usize, usize), usize>, // (chrom_index, site_pos1) -> count
    file: OutputFile,
}
impl Endpoints {

    /// Initialize endpoint extraction, opening a file prior to processing 
    /// in case IO fails.
    pub fn new(w: &mut Workflow) -> Endpoints {
        w.cfg.set_u8_env(&[MIN_MAPQ]);
        w.cfg.set_usize_env(&[CLIP_TOLERANCE]);
        w.cfg.set_usize_env(&[EXPECTING_ENDPOINT_RE_SITES]);
        let collecting = *w.cfg.get_bool(EXPECTING_ENDPOINT_RE_SITES);
        if collecting {
            w.ctrs.add_counters(&[
                (N_ENDS,      "read endpoints processed"),
                (N_HIGH_QUAL, format!(
                    "candidate endpoints with MAPQ >= {} and clip <= {}", 
                    w.cfg.get_u8(MIN_MAPQ), w.cfg.get_usize(CLIP_TOLERANCE)).as_str()
                ),
                (N_ENDPOINTS, "unique endpoint reference positions nominated as RE sites"),
            ]);
        }
        Endpoints {
            collecting,
            counter: HashMap::new(), // (chrom_index, site_pos1) -> count, incompatible with mdi::keyed_counters
            file: OutputFile::open_env(&mut w.cfg, OBSERVED_ENDPOINTS_FILE),
        }
    }

    /// Process the alignments for one read (of a pair).
    pub fn extract(
        read_n: usize, 
        alns: &Vec<&mut SamRecord>, 
        w: &mut Workflow,
        tool: &mut Tool,
        is_unmerged_pair: bool,
    ) -> Result<(), Box<dyn Error>> {
        if !tool.endpoints.collecting { return Ok(()); }

        // get ONT adapter trims
        let trims: Vec<usize> = if *w.cfg.get_bool(IS_ONT) {
            if let Some(tl) = alns[0].get_tag_value("tl") {
                tl.split(',').map(|x| x.parse().unwrap()).collect()
            } else {
                vec![0, 0]
            }
        } else {
            vec![0, 0]
        };

        // always attempt to process the 5' end of all reads
        Self::process_endpoint(alns[0], &trims, 5, w, tool)?;
        if read_n == 2 { return Ok(()); }

        // attempt to process read 3' ends if informative, i.e., inferred to be from end-to-end read
        // NOT necessarily on the same chromosome or part of the same contig as the 5' end
        if *w.cfg.get_bool(IS_END_TO_END_READ) ||    // PacBio or other platforms that guarantee end-to-end reads
            (
                *w.cfg.get_bool(IS_ONT) &&           // complete ONT reads, must be trimmed at 3' end to be confident end-to-end
                trims[TRIM3] > 0
            ) ||
            (
                *w.cfg.get_bool(IS_PAIRED_READS) &&  // merged paired reads
                !is_unmerged_pair
            ){
            Self::process_endpoint(alns[alns.len() - 1], &trims, 3, w, tool)?;
        }
        Ok(())
    }

    // check an informative endpoint for sufficient alignment quality
    fn process_endpoint(
        aln: &SamRecord, 
        trims: &Vec<usize>, 
        end: u8, 
        w: &mut Workflow,
        tool: &mut Tool,
    ) -> Result<(), Box<dyn Error>> {
        w.ctrs.increment(N_ENDS);

        // check for an informative endpoint with sufficient alignment quality
        // implicitly rejects unmapped with MAPQ=0
        if !tool.chroms.is_canonical(&aln.rname) || aln.mapq < *w.cfg.get_u8(MIN_MAPQ) { 
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
        if clip_len > *w.cfg.get_usize(CLIP_TOLERANCE) { 
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
        site_pos1 = (site_pos1 as i32 + clip_dir * clip_len as i32) as usize; // concern for possible overflow

        // use a typical trim value to adjust clips on ONT trim failures, to account for adapter bases still present
        // only applicable to 5' ends here, untrimmed 3' ends were filtered above as uninformative
        if *w.cfg.get_bool(IS_ONT) && end == 5 && trims[TRIM5] == 0 {
            site_pos1 = (site_pos1 as i32 - clip_dir * ONT_ADAPTER_LEN5) as usize;
        }

        // keep a tally of all genome positions nominated as RE sites
        *tool.endpoints.counter.entry((tool.chroms.index[&aln.rname], site_pos1)).or_insert(0) += 1;
        Ok(())
    }

    /// Report all unique observed endpoint positions with counts.
    pub fn write(
        w: &mut Workflow,
        tool: &mut Tool,
    ) -> Result<(), Box<dyn Error>> {
        if !tool.endpoints.collecting { return Ok(()); }
        w.log.print("printing endpoint tallies");
        let mut sorted_keys: Vec<(usize, usize)> = tool.endpoints.counter.keys().cloned().collect();
        sorted_keys.sort_unstable();
        for (chrom_index, site_pos1) in sorted_keys.iter() {
            w.ctrs.increment(N_ENDPOINTS);
            tool.endpoints.file.write_record(vec!(
                tool.chroms.rev_index[chrom_index].clone(),
                site_pos1.to_string(),
                tool.endpoints.counter[&(*chrom_index, *site_pos1)].to_string()
            ));
        }
        tool.endpoints.file.close();
        Ok(())
    }
}

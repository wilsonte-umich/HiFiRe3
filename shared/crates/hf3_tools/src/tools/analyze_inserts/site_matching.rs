//! Perform site matching for alignment endpoints.

// dependencies
use rustc_hash::FxHashMap;
use rayon::prelude::*;
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use mdi::InputFile;
use genomex::sam::{SamRecord, flag};
use crate::junctions::JxnFailureFlag;
use crate::sites::{SiteMatch, SiteMatches};

// constants
pub_key_constants!{
    // from environment variables
    ENZYME_NAME
    BLUNT_RE_TABLE
    OVERHANG5_RE_TABLE
    EXPECTING_ENDPOINT_RE_SITES
    REJECTING_JUNCTION_RE_SITES // can be true if EXPECTING_ENDPOINT_RE_SITES is false
    REJECT_JUNCTION_DISTANCE
    FILTERING_SITES_FILE 
    // counter keys
    N_SITES
    N_MATCHES_BY_OUTCOME
    // site match outcomes
    EXACT_SITE_MATCH
    INEXACT_SITE_MATCH
    TELOMERIC_POSITION
    NO_SITES_ON_CHROM
}

/// SiteMatcher matches alignment endpoints to RE filtering sites, i.e.,
/// site found by either in silico digestion of a genome or by examination
/// of read 5' endpoints.
pub struct SiteMatcher {
    pub is_matching_sites: bool, // false if not matching sites in ligFree or tagFree
    enzyme_name:           String,
    filtering_sites_file:  String,
    filtering_sites:       FxHashMap<String, Vec<u32>>, // chrom -> [site_pos1]
    n_sites:               usize,
    correction5:           u32, // 5' correction for site matching with 5' overhanging enzymes
    reject_junction_distance: i32,
}
impl SiteMatcher {
    /* ---------------------------------------------------------------------------
    initialize
    ---------------------------------------------------------------------------- */
    /// Initialize a new SiteMatcher.
    pub fn new(w: &mut Workflow) -> SiteMatcher {
        w.cfg.set_bool_env(&[EXPECTING_ENDPOINT_RE_SITES, REJECTING_JUNCTION_RE_SITES]);
        if !w.cfg.get_bool(REJECTING_JUNCTION_RE_SITES) {
            return SiteMatcher {
                is_matching_sites:      false,
                enzyme_name:            "".to_string(),
                filtering_sites_file:   "".to_string(),
                filtering_sites:        FxHashMap::default(),
                n_sites:                0,
                correction5:            0,
                reject_junction_distance: 0,
            };
        }
        w.cfg.set_string_env(&[ENZYME_NAME, BLUNT_RE_TABLE, OVERHANG5_RE_TABLE, 
                                     FILTERING_SITES_FILE]);
        w.cfg.set_u32_env( &[REJECT_JUNCTION_DISTANCE]);
        w.ctrs.add_counters(&[
            (N_SITES, &format!("{} filtering sites being matched against", w.cfg.get_string(ENZYME_NAME))),
        ]);
        w.ctrs.add_keyed_counters(&[
            (N_MATCHES_BY_OUTCOME, "RE site matching outcomes"),
        ]);
        let mut site_matcher = SiteMatcher{
            is_matching_sites:      true,
            enzyme_name:            w.cfg.get_string(ENZYME_NAME).to_string(),
            filtering_sites_file:    w.cfg.get_string(FILTERING_SITES_FILE).to_string(),
            filtering_sites:        FxHashMap::default(),
            n_sites:                0,
            correction5:            0,
            reject_junction_distance: *w.cfg.get_u32(REJECT_JUNCTION_DISTANCE) as i32,
        };
        site_matcher.load_filtering_sites();
        w.ctrs.add_to(N_SITES, site_matcher.n_sites);
        site_matcher.set_correction5(&w);
        site_matcher
    }

    // load the filtering sites from the assembled file
    fn load_filtering_sites(&mut self) {
        eprint!("loading RE filtering sites from {}", self.filtering_sites_file);
        for line in InputFile::get_lines(&self.filtering_sites_file, false) {
            let parts: Vec<&str> = line.split('\t').collect(); // chrom,sitePos1,inSilico,nObserved
            let chrom = parts[0];
            if chrom == "chrom" { continue; } // skip header line
            let site_pos1: u32 = parts[1].parse().unwrap_or(0);
            self.filtering_sites
                .entry(chrom.to_string())
                .or_insert_with(Vec::new)
                .push(site_pos1);
        }
        // ensure that each chrom's sites are sorted by site_pos1
        self.n_sites = 0;
        for sites in self.filtering_sites.values_mut() {
            sites.par_sort_unstable();
            self.n_sites += sites.len();
        }
    }
    // load the RE site metadata to properly handle blunt vs. 5' overhanging sites
    // enzyme  strand  cut_site regex   offset  CpG_priority
    // EcoRV   0       GATATC   GATATC  3       4  
    // NcoI    0       CCATGG   CCATGG  1       5  
    // thus:
    // blunt cutter: correction5 = 0
    //            *      sitePos1
    //        --3 5--    top strand alignment
    //  EcoRV GAT^ATC
    //        CTA^TAG
    //        --5 3--    bottom strand alignment
    // 5' overhanging cutter: correction5 = siteLength - 2 * offset
    //          *        sitePos1
    //      --3 5--      top strand alignment
    //  NcoI  C^CATG G   offset = 1, correction5 = 6 - 2 * 1 = 4
    //        G GTAC^C
    //           --5 3-- bottom strand alignment
    // correction5 is applied to alignment pos1 on reverse strand by apply_filter.pl when querying for closest site
    fn set_correction5(&mut self, w: &Workflow) {
        let mut correction5: Option<u32> = None;
        for table in &[w.cfg.get_string(BLUNT_RE_TABLE), w.cfg.get_string(OVERHANG5_RE_TABLE)] {
            // enzyme,strand,cut_site,regex,offset,CpG_priority,high_fidelity,site_length
            for line in InputFile::get_lines(table, true) {
                let parts: Vec<&str> = line.split(',').collect();
                let enzyme = parts[0].trim();
                if enzyme == self.enzyme_name {
                    if table == &w.cfg.get_string(BLUNT_RE_TABLE) {
                        correction5 = Some(0);
                    } else {
                        let offset: u32 = parts[4].parse().unwrap();
                        let site_length: u32 = parts[7].parse().unwrap();
                        correction5 = Some(site_length - 2 * offset); // account for staggered cleaved bonds with 5' overhanging enzymes
                    }
                    break;
                }
            }
        }
        self.correction5 = if let Some(correction5) = correction5 {
            correction5
        } else {
            panic!("could not find enzyme {} in RE tables", self.enzyme_name);
        };
    }
    /* ---------------------------------------------------------------------------
    genome position to RE site matching
    ---------------------------------------------------------------------------- */
    /// Match a specific position on a chromosome to the closest RE filtering site.
    pub fn find_closest_site(
        &self, 
        chrom:    &str, 
        qry_pos1: u32,
        w:        &mut Workflow, 
    ) -> SiteMatch {
        if let Some(sites) = self.filtering_sites.get(chrom) {
            match sites.binary_search(&qry_pos1) {
                Ok(site_index0) => { // exact match
                    w.ctrs.increment_keyed(N_MATCHES_BY_OUTCOME, EXACT_SITE_MATCH);
                    SiteMatch {
                        index1: site_index0 as i32 + 1, // 1-based index
                        pos1: sites[site_index0],
                        distance: 0,
                    }
                },
                Err(next_site_index0) => { // no exact match
                    // return the closest site with index1 and distance signed negative if qry_pos1 < site.pos1
                    if next_site_index0 > 0 && next_site_index0 < sites.len() {
                        w.ctrs.increment_keyed(N_MATCHES_BY_OUTCOME, INEXACT_SITE_MATCH);
                        let dist_low = qry_pos1 - sites[next_site_index0 - 1];
                        let dist_high = sites[next_site_index0] - qry_pos1;
                        if dist_low < dist_high {
                            SiteMatch {
                                index1: next_site_index0 as i32, // 1-based index, i.e., index1 == next_site_index0 - 1 + 1
                                pos1: sites[next_site_index0 - 1],
                                distance: dist_low as i32,
                            }
                        } else {
                            SiteMatch {
                                index1: -(next_site_index0 as i32 + 1), // 1-based index
                                pos1: sites[next_site_index0],
                                distance: -(dist_high as i32),
                            }
                        }
                    // telomeric sites on a chromosome are considered unusable and return a null site
                    } else {
                        w.ctrs.increment_keyed(N_MATCHES_BY_OUTCOME, TELOMERIC_POSITION);
                        SiteMatch::new()
                    }
                },
            }
        } else {
            w.ctrs.increment_keyed(N_MATCHES_BY_OUTCOME, NO_SITES_ON_CHROM);
            SiteMatch::new()
        }
    }
    /// Return the two `SiteMatch`es for SamRecord alignment as (5' site, 3' site)
    /// in read order.
    pub fn match_aln_sites(
        &self, 
        aln: &SamRecord,
        w:   &mut Workflow, 
    ) -> (SiteMatch, SiteMatch) {
        if self.is_matching_sites {
            if aln.check_flag_any(flag::REVERSE){
                (
                    self.find_closest_site(&aln.rname, aln.get_end1() + 1 - self.correction5, w),
                    self.find_closest_site(&aln.rname, aln.pos1 - self.correction5, w),
                )
            } else {
                (
                    self.find_closest_site(&aln.rname, aln.pos1, w),
                    self.find_closest_site(&aln.rname, aln.get_end1() + 1, w),
                )
            }
        } else {
            (SiteMatch::new(), SiteMatch::new())
        }
    }
    /* ---------------------------------------------------------------------------
    project alignment 3' end to next RE site
    ---------------------------------------------------------------------------- */
    /// Return the projected 3' SiteMatch for a SamRecord alignment.
    pub fn get_projection(
        &self, 
        aln:       &SamRecord,
        mut site3: SiteMatch, // provided as the actual aln 3' site to be projected
    ) -> SiteMatch {
        // find the next site in the direction of the alignment strand
        // actual and projected site indices could be the same
        // no sign on index, pretend an exact match
        //          +          -
        //  ----|------------------|-----
        //      |------------->~~~~~      first three DO NOT need to be projected, closest site is the projected site
        //      |-------------------> (small overrun)
        //      ~~~~<--------------|
        //      |------>~~~~~~~~~~~~      last three DO need to be projected, closet site is one away from the projected site
        //      ~~~~~~~~~~~~~<-----|       
        //     <-------------------|
        let chrom_sites = self.filtering_sites.get(&aln.rname).unwrap();
        if aln.check_flag_any(flag::REVERSE){
            if site3.index1 < 0 && 
               site3.distance.abs() > self.reject_junction_distance && // more permissive than ACCEPT_ENDPOINT_DISTANCE
               site3.index1.abs() > 1 {
                site3.index1 = site3.index1.abs() - 1;
            }
        } else {
            let n_sites_on_chrom = chrom_sites.len() as i32;
            if site3.index1 > 0 && 
               site3.distance.abs() > self.reject_junction_distance &&
               site3.index1 < n_sites_on_chrom {
                site3.index1 += 1; 
            }
        }
        site3.index1 = site3.index1.abs(); // projected site has no sign, implies exact match
        site3.pos1 = chrom_sites[(site3.index1 - 1) as usize];
        site3.distance = 0;
        site3
    }
    /* ---------------------------------------------------------------------------
    perform chimeric junction assessment based on RE site matches
    ---------------------------------------------------------------------------- */
    /// Assess whether a junction has either end within REJECT_JUNCTION_DISTANCE
    /// of an RE filtering site, returning JxnFailureFlag::SiteMatch if the junction 
    /// should be rejected.
    pub fn check_jxn(
        &self,  
        aln_sites: &[SiteMatches],
        aln5_i:    usize, // the alignment 5' to the junction on the read (its 3' end is the breakpoint)
        aln3_i:    usize,
        w: &mut Workflow,
    ) -> JxnFailureFlag {
        if !self.is_matching_sites { return JxnFailureFlag::None; }
        let dist_prox = &aln_sites[aln5_i].site3.distance;
        let dist_dist = &aln_sites[aln3_i].site5.distance;
        // either breakpoint node can fail the junction site match test
        let failed = dist_prox.abs() <= self.reject_junction_distance ||
                           dist_dist.abs() <= self.reject_junction_distance;
        if failed {
            w.ctrs.increment_keyed(super::N_JXNS_BY_REASON, super::JXN_FAIL_SITE_MATCH);
            JxnFailureFlag::SiteMatch
        } else {
            JxnFailureFlag::None
        }
    }
}

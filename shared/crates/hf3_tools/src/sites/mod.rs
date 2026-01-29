//! Support site matching for alignment endpoints.

// dependencies
use rust_htslib::bam::record::Record as BamRecord;
use genomex::bam::tags;
use crate::formats::hf3_tags::*;

/// SiteMatch structure to declare the best match of a genome 
/// position/strand to an RE filtering site.
#[derive(Clone)]
pub struct SiteMatch {
    // signed 1-based index in filtering_sites_file of the closest RE site to query.pos1
    // sign = direction from query.pos1 to site.pos1, + means query.pos1 >= site.pos1
    // 0 = no site matched
    // indices reset to 1 at the beginning of each chromosome
    pub index1:   i32, 
    // 1-based site coordinate on the chromosome
    pub pos1:     u32, 
    // query.pos1 - site.pos1
    pub distance: i32, 
}
impl SiteMatch {
    /// Return a new empty SiteMatch.
    pub fn new() -> SiteMatch {
        SiteMatch {
            index1:   0,
            pos1:     0,
            distance: 0,
        }
    }
}

/// SiteMatches structure to hold 5', 3', and projected 3' site matches.
pub struct SiteMatches {
    pub site5: SiteMatch,
    pub site3: SiteMatch,
    pub proj3: SiteMatch,
}
impl SiteMatches {
    /* ---------------------------------------------------------------------------
    pack and unpack sam tags
    ---------------------------------------------------------------------------- */
    /// Pack site matches into a string representation of i-type closest sites
    /// as (site5, site3, proj3) for use as a value for tag prefix 'sc:B:i,'.
    pub fn to_tag(&self) -> String {
        [   
            self.site5.index1.to_string(), // 5' site
            self.site5.pos1.to_string(),
            self.site5.distance.to_string(),
            self.site3.index1.to_string(), // 3' site
            self.site3.pos1.to_string(),
            self.site3.distance.to_string(),
            self.proj3.index1.to_string(), // projected 3' site
            self.proj3.pos1.to_string(),
            self.proj3.distance.to_string(),
        ].join(",")
    }
    
    /// Unpack site matches from an alignment's 'sc:B:i,' tag 
    /// into SiteMatches{site5, site3, proj3}.
    pub fn from_bam_record(aln: &BamRecord) -> SiteMatches {
        let vals = tags::get_tag_i32_vec_opt(aln, CLOSEST_SITES);
        if let Some(vals) = vals {
            let site5 = SiteMatch {
                index1:   vals[0] as i32,
                pos1:     vals[1] as u32,
                distance: vals[2] as i32,
            };
            let site3 = SiteMatch {
                index1:   vals[3] as i32,
                pos1:     vals[4] as u32,
                distance: vals[5] as i32,
            };
            let proj3 = SiteMatch {
                index1:   vals[6] as i32,
                pos1:     vals[7] as u32,
                distance: vals[8] as i32,
            };
            SiteMatches{site5, site3, proj3}
        } else {
            SiteMatches {
                site5: SiteMatch::new(),
                site3: SiteMatch::new(),
                proj3: SiteMatch::new(),
            }
        }
    }
}

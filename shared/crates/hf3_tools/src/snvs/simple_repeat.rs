/// Support for masking SNV/indels in simple repeats, which are error prone
/// in PacBio sequencing. 

// imports
use std::cmp::Ordering;
use rustc_hash::FxHashMap;
use serde::Deserialize;
use mdi::InputCsv;
use crate::snvs::*;

// constants
const MIN_HOMOPOLYMER_LENGTH: u32 = 6; // TODO: expose as option?
const SEARCH_PADDING: u32 = 1;

/// Structure for collecting simple repeat spans. All variants in these spans
/// are disallowed during rare SNV/indel calling.
#[derive(PartialEq, Eq, PartialOrd, Ord)]
pub struct SimpleRepeat {
    pub start0: u32,  // 0-based, inclusive
    pub end1:   u32,  // 0-based, exclusive (open interval)
}

/// Structure for loading RepeatMasker simple repeats from BED.
#[derive(Deserialize)]
struct RmskRecord {
    pub chrom:  String,
    pub start0: u32,  // 0-based, inclusive
    pub end1:   u32,  // 0-based, exclusive (open interval)
    pub _unit:  String,
}

/// Structure for loading TandemRepeatFinder simple repeats from BED.
#[derive(Deserialize)]
struct TrfRecord {
    pub chrom:  String,
    pub start0: u32,  // 0-based, inclusive
    pub end1:   u32,  // 0-based, exclusive (open interval)
    pub _period:      String,
    pub _copy_number: String,
}

/// SimpleRepeats is a searchable collection of simple repeat spans collapsed
/// from three methods: HiFiRe3 homopolymer run finding, RepeatMasker
/// simple_repeat annotation, and TandemRepeatFinder.
pub struct SimpleRepeats{
    repeats:     Vec<SimpleRepeat>,
    by_fragment: FxHashMap<ReFragment, ChromLength>,
}
impl SimpleRepeats {

    /// Fill the simple repeat spans on a specific chromosome.
    pub fn new(
        tool:  &SnvAnalysisTool,
        chrom: &str,
        rmsk_simple_repeats_bed: &str,
        trf_simple_repeats_bed:  &str,
    ) -> Self {
        let mut tmp_repeats: Vec<SimpleRepeat> = Vec::with_capacity(1_000_000);

        // load homopolymer runs; this is more sensitive than RepeatMasker
        let tid = tool.fa.fai().tid(chrom).expect("Cannot find homopolymer {chrom}");
        let chrom_view = tool.fa.view_tid(tid).expect("Cannot find homopolymer {chrom}");
        let mut base: Option<u8> = None;
        let mut start0 = 0_u32;
        let mut len = 0_u32;
        for (pos0, base_ref) in chrom_view.bases().enumerate() {
            let qry_base = base_ref.to_ascii_uppercase();
            match base {
                Some(base) if qry_base == base => {
                    len += 1;
                },
                _ => {
                    if len >= MIN_HOMOPOLYMER_LENGTH {
                        tmp_repeats.push(
                            SimpleRepeat { start0, end1: pos0 as u32 }
                        );
                    }
                    base = Some(qry_base);
                    start0 = pos0 as u32;
                    len = 1;
                },
            }
        }
        // ignore the last run as inconsequential

        // further load simple repeats called by RepeatMasker ...      
        let mut bed= InputCsv::open_file(rmsk_simple_repeats_bed, b'\t', false);
        for result in bed.deserialize::<RmskRecord>() {
            let record = result.expect("Could not extract record from RMSK row.");
            if record.chrom == chrom {
                tmp_repeats.push(
                    SimpleRepeat { start0: record.start0, end1: record.end1 }
                );                
            }
        }

        // ... and also by Tandem Repeat finder 
        // !! TODO: at present, this input file must be created manually from UCSC !!
        let path = std::path::Path::new(trf_simple_repeats_bed);
        if path.exists(){
            let mut bed= InputCsv::open_file(trf_simple_repeats_bed, b'\t', false);
            for result in bed.deserialize::<TrfRecord>() {
                let record = result.expect("Could not extract record from TRF row.");
                if record.chrom == chrom {
                    tmp_repeats.push(
                        SimpleRepeat { start0: record.start0, end1: record.end1 }
                    );                
                }
            }
        }

        // sort and collapse the two types of entries together for binary search
        tmp_repeats.sort_unstable();
        let mut simple_repeats: Vec<SimpleRepeat> = Vec::new();
        for tmp_repeat in tmp_repeats {
            match simple_repeats.last_mut() {
                Some(last) if tmp_repeat.start0 <= last.end1 => {
                    last.end1 = last.end1.max(tmp_repeat.end1);
                }
                _ => simple_repeats.push(tmp_repeat),
            }
        }
        let mut value = Self {
            repeats:     simple_repeats,
            by_fragment: FxHashMap::default(),
        };
        value.by_fragment.reserve(4096);
        value
    }

    /// Use a binary search to determine if a variant span overlaps a simple 
    /// repeat.
    /// 
    /// Use a 1-bp padding to reject immediately adjacent variants also.
    pub fn binary_search(
        &self, 
        mut start0: ChromPos0, 
        len: ChromLength
    ) -> bool {
        let end1 = start0 + len + SEARCH_PADDING;
        start0 = start0.saturating_sub(SEARCH_PADDING);
        let result = self.repeats.binary_search_by(|simple_repeat| {
            if simple_repeat.end1 <= start0 {
                Ordering::Less
            } else if simple_repeat.start0 >= end1 {
                Ordering::Greater
            } else {
                Ordering::Equal
            }
        });
        result.is_ok()
    }

    /// Return the number of bases in a genome span that are masked as simple
    /// repeats.
    /// 
    /// Use a 1-bp padding to reject immediately adjacent variants also. 
    pub fn get_n_repeat_bases(
        &mut self, 
        re_fragment: &ReFragment, 
    ) -> ChromLength {
        if let Some(n_repeat_bases) = self.by_fragment.get(re_fragment) {
            *n_repeat_bases
        } else {
            let n_repeat_bases = self.get_base_count(re_fragment);
            self.by_fragment.insert(*re_fragment, n_repeat_bases);
            n_repeat_bases
        }
    }
    fn get_base_count(
        &self, 
        re_fragment: &ReFragment, 
    ) -> ChromLength {

        // get the first repeat within or after the fragment
        let result = self.repeats.binary_search_by(|simple_repeat| {
            if simple_repeat.end1 <= re_fragment.start0 {
                Ordering::Less
            } else if simple_repeat.start0 > re_fragment.start0 {
                Ordering::Greater
            } else { // rare case where repeat starts at the beginning of the RE fragment
                Ordering::Equal 
            }
        });

        // scan from the leftmost found repeat span
        match result {
            Err(i0) => {
                self.scan_repeats(re_fragment, i0)
            },
            Ok(i0) => {
                self.scan_repeats(re_fragment, i0)
            },
        }
    }
    fn scan_repeats(
        &self, 
        re_fragment: &ReFragment, 
        mut i0: usize,
    ) -> ChromLength {
        let mut n_repeat_bases: ChromLength = 0;
        let max_i0 = self.repeats.len() - 1;
        while i0 <= max_i0 && self.repeats[i0].start0 < re_fragment.end1 {
            let overlap_start0 = re_fragment.start0.max(self.repeats[i0].start0);
            let overlap_end1 = re_fragment.end1.min(self.repeats[i0].end1);
            n_repeat_bases += overlap_end1.saturating_sub(overlap_start0) + 2;
            i0 += 1;
        }
        n_repeat_bases
    }
}

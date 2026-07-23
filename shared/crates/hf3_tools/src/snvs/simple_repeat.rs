/// Support for masking SNV/indels in simple repeats, which are error prone
/// in PacBio sequencing. 

// dependencies
use std::cmp::Ordering;
use serde::Deserialize;
use mdi::InputCsv;
use crate::snvs::SnvAnalysisTool;

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
pub struct SimpleRepeats(Vec<SimpleRepeat>);
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
        for simple_repeat in tmp_repeats {
            match simple_repeats.last_mut() {
                Some(last) if simple_repeat.start0 <= last.end1 => {
                    last.end1 = last.end1.max(simple_repeat.end1);
                }
                _ => simple_repeats.push(simple_repeat),
            }
        }
        Self(simple_repeats)
    }

    /// Use a binary search to determine if a variant span overlaps a simple 
    /// repeat.
    /// 
    /// Use a 1-bp padding to reject immediately adjacent variants also.
    pub fn binary_search(&self, mut start0: u32, len: u32) -> bool {
        let end1 = start0 + len + SEARCH_PADDING;
        start0 = start0.saturating_sub(SEARCH_PADDING);
        let result = self.0.binary_search_by(|simple_repeat| {
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
}

//! Support for creating encoded fragment-level files for downstream use.

// imports
use rustc_hash::FxHashMap;
use mdi::OutputCsv;
use crate::snvs::*;

// constants
const MATCH:   char = '='; // homoduplex reference base matchs
const MASKED:  char = 'N'; // heteroduplex bases considered untrustworthy for variant calling
const ALT_DEL: char = '-'; // fully recorded in the reference-justified encoded "read"
const ALT_DEL_MASKED: char = '?';
const ALT_INS: char = '+'; // recorded in a separate encoding from match/subs/del
const ALT_INS_MASKED: char = '?'; 

/// AlignmentEncoding progressively builds the encoded bases of a single 
/// sequence (either a read or a reference sequence) aligned to its reference
/// (either the reference sequence or a haplotype consensus).
#[derive(Clone)]
pub struct AlignmentEncoding {
    re_fragment: ReFragment,
    n_match:     u32,
    n_alt:       u32,
    n_del:       u32,
    n_ins:       u32,
    n_masked:    u32,
    sample_bit:  SampleBit, // sample bit for this variant
    read_start0: SeqPos0, // target base where encoding starts, may be less than re_fragment.start0
    encoding:    String,
    insertion:   Vec<String>, // so insertions can be tracked on top of any previous base
    qname:       QName,
}
impl AlignmentEncoding {

    /// Create a new AlignmentEncoding object. One AlignmentEncoding is created
    /// per SnvChromWorker that is reset as needed per read.
    pub fn new() -> Self {
        AlignmentEncoding { 
            re_fragment: ReFragment { start0: 0, end1: 0 }, 
            n_match:     0,
            n_alt:       0,
            n_del:       0,
            n_ins:       0, // not included as part of n_bases
            n_masked:    0,
            sample_bit:  0,
            read_start0: 0,
            encoding:    String::with_capacity(256),
            insertion:   Vec::with_capacity(256), // each encodes the same number of ref bases
            qname:       String::with_capacity(64),
        }
    }

    /// Prepare an allocated AlignmentEncoding to build a new read encoding
    /// based on a reference chromosome alignment.
    pub fn prepare_read_on_ref(
        &mut self,
        re_fragment: &ReFragment,
        read: &ReadInstance,
    ) {
        self.re_fragment = *re_fragment;
        self.n_match  = 0;
        self.n_alt    = 0;
        self.n_del    = 0;
        self.n_ins    = 0;
        self.n_masked = 0;
        self.sample_bit  = read.sample_bit;
        self.read_start0 = read.ref_pos0.min(re_fragment.start0);
        self.encoding.clear();
        self.insertion.clear();
        self.qname.clear();
        self.qname.push_str(&read.qname); 
        if read.ref_pos0 > re_fragment.start0 {
            let len = (read.ref_pos0 - re_fragment.start0) as usize;
            self.encoding.extend(vec![MASKED; len]);
            self.insertion.push(format!("{}{}", MATCH, len));
        }
    }

    /// Prepare an allocated AlignmentEncoding to build a new read encoding
    /// based on a haplotype consensus alignment.
    pub fn prepare_read_on_hap(
        &mut self,
        re_fragment: &ReFragment,
        read: &ReadInstance,
        tgt_start0: usize,
    ) {
        self.re_fragment = *re_fragment;
        self.n_match  = 0;
        self.n_alt    = 0;
        self.n_del    = 0;
        self.n_ins    = 0;
        self.n_masked = 0;
        self.sample_bit  = read.sample_bit;
        self.read_start0 = 0; // unlike read_on_ref, read_on_hap cannot align to bases before the consensus starts
        self.encoding.clear();
        self.insertion.clear();
        self.qname.clear();
        self.qname.push_str(&read.qname); 
        if tgt_start0 > 0 {
            self.encoding.extend(vec![MASKED; tgt_start0]);
            self.insertion.push(format!("{}{}", MATCH, tgt_start0));
        }
    }
    
    /// Add a base span identical to target to a growing read encoding.
    pub fn add_identity(&mut self, len: u32) {
        self.n_match += len;
        let op = format!("{}{}", MATCH, len);
        self.encoding.push_str(&op);
        self.insertion.push(op);
    }
    /// Add a base substitution to a growing read encoding.
    pub fn add_substitution(&mut self, alt_base: &str, allowed: bool, low_qual: bool) {
        if allowed {
            self.n_alt += 1;
            if low_qual { 
                self.encoding.push_str(&alt_base.to_ascii_lowercase());
            } else {
                self.encoding.push_str(alt_base);
            }
        } else {
            self.n_masked += 1;
            self.encoding.push(MASKED);
        }
        self.insertion.push(format!("{}{}", MATCH, 1));
    }
    /// Add an insertion to a growing read encoding on the last known base.
    pub fn add_insertion(&mut self, allowed: bool, low_qual: bool) {
        self.n_ins += 1; // number of insertion events, not bases
        if let Some(prev) = self.insertion.pop() { // for safety; never expection ins as first op
            let op = prev.chars().next().unwrap(); // cannot be empty
            if op == MATCH { // for safety; never expect two insertions in a row
                let len = &mut prev[1..].parse::<usize>().unwrap();
                *len -= 1; // reduce the previous non-insertion run by one
                if *len > 0 {
                    self.insertion.push(format!("{}{}", MATCH, len));
                }   
                self.insertion.push(if allowed && !low_qual { // than add this insertion on that one base
                    ALT_INS.to_string()
                } else {
                    ALT_INS_MASKED.to_string()
                });
            }
        }
        // do nothing to self.encoding
    }
    /// Add a deletion to a growing read encoding.
    pub fn add_deletion(&mut self, len: u32, allowed: bool, low_qual: bool) {
        let op = if allowed && !low_qual  {
            self.n_del += len;
            ALT_DEL
        } else {
            self.n_masked += len;
            ALT_DEL_MASKED
        };
        self.encoding.extend(vec![op; len as usize]);
        self.insertion.push(format!("{}{}", MATCH, len));
    }
}

/// EncodingMetadata reports summary results of read encoding accumulated
/// over all reads.
pub struct EncodingMetadata {
    pub n_unique_spans:  usize,
    pub n_reads:         usize,
    pub n_ref_bases:     usize,
    pub n_read_bases:    usize,
    pub n_match:         usize,
    pub n_alt:           usize,
    pub n_del:           usize,
    pub n_ins:           usize,
    pub n_masked:        usize,
}
impl EncodingMetadata {
    fn new(n_unique_spans: usize) -> Self {
        EncodingMetadata {
            n_unique_spans,
            n_reads:         0,
            n_ref_bases:     0,
            n_read_bases:    0,
            n_match:         0,
            n_alt:           0,
            n_del:           0,
            n_ins:           0,
            n_masked:        0,
        }
    }
}

/// A fully assembled encoding over multiple aggregated reads as written to file  
/// for downstream use. 
#[derive(Serialize)]
struct EncodingRecord {
    chrom_index:  ChromIndex1,
    re_fragment:  ReFragment,
    #[serde(serialize_with = "serialize_haplotype")]
    haplotype:    Haplotype,
    n_reads:      usize,
    n_match:      u32, // aggregated valued over all reads
    n_alt:        u32,
    n_del:        u32,
    n_ins:        u32,
    n_masked:     u32,
    sample_bits:  SampleBits,
    n_matches:    CommaDelimited, // comma-delimited values from the invidual reads
    n_alts:       CommaDelimited,
    n_dels:       CommaDelimited,
    n_inss:       CommaDelimited,
    n_maskeds:    CommaDelimited,
    sample_bitss: CommaDelimited,
    read_start0s: CommaDelimited,
    encodings:    CommaDelimited,
    insertions:   CommaDelimited,
    qnames:       CommaDelimited,
    seq:          UppercaseACGTN, // of the top ref strand from start0 to end1
    hap_vs_ref:   String, // encoding of haplotype consensus differences relative to reference
}

/// AlignmentEncodings aggregates AlignmentEncodings per ReFragment.
pub struct AlignmentEncodings{
    encodings: FxHashMap<(ReFragment, Haplotype), Vec<AlignmentEncoding>>
}
impl AlignmentEncodings {

    /// Create a new empty AlignmentEncodings object. On encoding tally is 
    /// instantiated per SnvChromWorker.
    pub fn new() -> Self{
        let mut encodings = FxHashMap::default();
        encodings.reserve(4096); // about 16 Mb of ReFragments to start out
        Self{ encodings }
    }

    /// Add a completed alignment encoding to the map.
    pub fn insert(
        &mut self, 
        re_fragment: &ReFragment, 
        haplotype:   Haplotype,
        encoding:    AlignmentEncoding
    ) {
        self.encodings.entry((*re_fragment, haplotype))
            .or_insert_with(|| Vec::with_capacity(8))
            .push(encoding);          
    }

    /// Sort and write a set of AlignmentEncodings to a temporary file for the
    /// working chromosome.
    pub fn write_sorted(
        tool:      &SnvAnalysisTool,
        worker:    &SnvChromWorker,
        haplotype_consensuses: &mut HaplotypeConsensuses,
        encodings: &AlignmentEncodings, // since there are multiple encodings files
        encodings_file_path: &str,
    ) -> EncodingMetadata {
        let mut csv = OutputCsv::open_csv(
            encodings_file_path, 
            b'\t', 
            false, 
            Some(tool.n_cpu),
        );
        let mut read_sets = encodings.encodings.keys()
            .filter_map(|rs|{
                let excluded  =  tool.exclusions.pos_in_region(&worker.chrom, rs.0.start0 + 1);
                let on_target = !tool.targets.has_data || 
                                       tool.targets.pos_in_region(&worker.chrom, rs.0.start0 + 1);
                if !excluded && on_target { Some(rs) } else { None }
            }).collect::<Vec<_>>();
        read_sets.sort_unstable();
        let mut md = EncodingMetadata::new(read_sets.len());
        for rs in read_sets {
            let encodings = &encodings.encodings[rs];
            let n_reads = encodings.len();
            let (hap_seq, hap_vs_ref) = haplotype_consensuses.cache
                .get_mut(rs)
                .expect("Failed to get haplotype seq from cache.");
            let record = EncodingRecord {
                chrom_index: worker.chrom_index,
                re_fragment: rs.0,
                haplotype:   rs.1,
                n_reads:     n_reads,

                n_match:     encodings.iter().map(|r| r.n_match).sum::<u32>(),
                n_alt:       encodings.iter().map(|r| r.n_alt).sum::<u32>(),
                n_del:       encodings.iter().map(|r| r.n_del).sum::<u32>(),
                n_ins:       encodings.iter().map(|r| r.n_ins).sum::<u32>(),
                n_masked:    encodings.iter().map(|r| r.n_masked).sum::<u32>(),

                sample_bits: encodings.iter().fold(0_u32, |acc, rd| acc | rd.sample_bit),

                n_matches:   encodings.iter().map(|r| r.n_match.to_string())
                                .collect::<Vec<String>>().join(","),
                n_alts:      encodings.iter().map(|r| r.n_alt.to_string())
                                .collect::<Vec<String>>().join(","),
                n_dels:      encodings.iter().map(|r| r.n_del.to_string())
                                .collect::<Vec<String>>().join(","),
                n_inss:      encodings.iter().map(|r| r.n_ins.to_string())
                                .collect::<Vec<String>>().join(","),
                n_maskeds:   encodings.iter().map(|r| r.n_masked.to_string())
                                .collect::<Vec<String>>().join(","),

                sample_bitss: encodings.iter().map(|r| r.sample_bit.to_string())
                                .collect::<Vec<String>>().join(","),
                read_start0s: encodings.iter().map(|r| r.read_start0.to_string())
                                .collect::<Vec<String>>().join(","),                   

                encodings:    encodings.iter().map(|r| r.encoding.clone())
                                .collect::<Vec<String>>().join(","),
                insertions:   encodings.iter().map(|r| r.insertion.join(""))
                                .collect::<Vec<String>>().join(","),

                qnames:       encodings.iter().map(|r| r.qname.clone())
                                .collect::<Vec<String>>().join(","),

                seq:        std::mem::take(hap_seq),
                hap_vs_ref: std::mem::take(hap_vs_ref).unwrap_or_else(|| "NA".to_string()),

            };
            csv.serialize(&record);

            let fragment_len = (rs.0.end1 - rs.0.start0) as usize;
            md.n_reads      += n_reads;
            md.n_ref_bases  += fragment_len;
            md.n_read_bases += n_reads * fragment_len;
            md.n_match      += record.n_match as usize;
            md.n_alt        += record.n_alt as usize;
            md.n_del        += record.n_del as usize;
            md.n_ins        += record.n_ins as usize;
            md.n_masked     += record.n_masked as usize;
        }
        md
    }
}

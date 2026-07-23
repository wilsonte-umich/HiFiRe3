//! Support for creating an encoded read-level file for visualization in the app.

// // dependencies
// use std::str::from_utf8_unchecked;
// use rustc_hash::FxHashMap;
// use rust_htslib::bam::Record as BamRecord;
// use genomex::sequence::rc_acgtn_str;
// use genomex::bam::tags as bam_tags;
// use crate::formats::hf3_tags::*;
// use super::*;

// // constants
// const MATCH:   char = '='; // homoduplex reference base matchs
// const MASKED:  char = 'N'; // heteroduplex bases considered untrustworthy for variant calling
// const ALT_DEL: char = '-'; // fully recorded in the reference-justified encoded "read"
// const ALT_DEL_MASKED: char = '?';
// const ALT_INS: char = '+'; // recorded in a separate encoding from match/subs/del
// const ALT_INS_MASKED: char = '?'; 

// /// A ReadEncoding progressively builds the encoded bases of a single read.
// pub struct ReadEncoding {
//     start0:   u32, // BED half-open coordinates on a known chromosome.
//     end1:     u32, // these are the fixed site positions the read matched to
//     n_match:  u32,
//     n_alt:    u32,
//     n_masked: u32,
//     n_del:    u32,
//     n_ins:    u32,
//     read_start0: u32, // ref base where encoding starts, may be less than start0
//     strand:      u8,
//     encoding:    String,
//     insertion:   Vec<String>, // so insertions can be tracked on top of any previous base
//     sample_bit:  u32,   // sample bits for this variant
//     qname:       String,
// }
// impl ReadEncoding {
//     /// Create a new ReadEncoding object.
//     pub fn new(
//         sample_bit: u32,
//         qname:    &[u8], 
//         ref_pos0: u32, 
//         mut site5_pos1: u32, 
//         mut site3_pos1: u32,
//         strand:   u8,
//     ) -> Self {
//         if strand == 1 { // reverse endpoints for rare reverse strand reads
//             (site5_pos1, site3_pos1) = (site3_pos1, site5_pos1);
//         }
//         let start0 = site5_pos1 - 1;
//         let mut encoding = ReadEncoding { 
//             start0, // encoding may start before start0, at read_start0
//             end1:      site3_pos1 - 1, 
//             n_match:   0,
//             n_alt:     0,
//             n_masked:  0,
//             n_del:     0,
//             n_ins:     0, // not included as part of n_bases
//             read_start0: ref_pos0.min(start0),
//             strand,
//             encoding:  String::with_capacity(200),
//             insertion: Vec::with_capacity(200), // each encodes the same number of ref bases
//             sample_bit,
//             qname:     unsafe { from_utf8_unchecked(qname).to_string() },
//         };
//         if ref_pos0 > start0 {
//             let len = (ref_pos0 - start0) as usize;
//             encoding.encoding.extend(vec![MASKED; len]);
//             encoding.insertion.push(format!("{}{}", MATCH, len));
//         }
//         encoding
//     }
//     /// Add a reference base span to a growing read encoding.
//     pub fn add_ref(&mut self, len: u32) {
//         self.n_match += len;
//         let op = format!("{}{}", MATCH, len);
//         self.encoding.push_str(&op);
//         self.insertion.push(op);
//     }
//     /// Add a base substitution to a growing read encoding.
//     pub fn add_subs(&mut self, alt_base: char, allowed: bool) {
//         self.encoding.push(if allowed {
//             self.n_alt += 1;
//             alt_base
//         } else {
//             self.n_masked += 1;
//             MASKED
//         });
//         self.insertion.push(format!("{}{}", MATCH, 1));
//     }
//     /// Add an insertion to a growing read encoding on the last known base.
//     pub fn add_ins(&mut self, allowed: bool) {
//         self.n_ins += 1; // number of insertion events, not bases
//         if let Some(prev) = self.insertion.pop() { // for safety; never expection ins as first op
//             let op = prev.chars().next().unwrap(); // cannot be empty
//             if op == MATCH { // for safety; never expect two insertions in a row
//                 let len = &mut prev[1..].parse::<usize>().unwrap();
//                 *len -= 1; // reduce the previous non-insertion run by one
//                 if *len > 0 {
//                     self.insertion.push(format!("{}{}", MATCH, len));
//                 }   
//                 self.insertion.push(if allowed { // than add this insertion on that one base
//                     ALT_INS.to_string()
//                 } else {
//                     ALT_INS_MASKED.to_string()
//                 });
//             }
//         }
//         // do nothing to self.encoding
//     }
//     /// Add a deletion to a growing read encoding.
//     pub fn add_del(&mut self, len: u32, allowed: bool) {
//         let op = if allowed {
//             self.n_del += len;
//             ALT_DEL
//         } else {
//             self.n_masked += len;
//             ALT_DEL_MASKED
//         };
//         self.encoding.extend(vec![op; len as usize]);
//         self.insertion.push(format!("{}{}", MATCH, len));
//     }
// }

// /// A fully assembled encoding over multiple aggregated reads to be serialized
// /// into the BED output.
// #[derive(Serialize)]
// struct EncodingRecord {
//     chrom_index: u8,
//     start0:    u32, // BED half-open coordinates
//     end1:      u32,
//     n_reads:   u32,
//     n_match:   u32, // summed valued over all reads
//     n_alt:     u32,
//     n_masked:  u32,
//     n_del:     u32,
//     n_ins:     u32,
//     n_reverse: u32,
//     n_matches: String, // comma-delimited values from the invidual reads
//     n_alts:    String,
//     n_maskeds: String,
//     n_dels:    String,
//     n_inss:    String,
//     read_start0s: String,
//     strands:      String,
//     encodings:    String,
//     insertions:   String,
//     sample_bits:  u32,
//     sample_bitss: String,
//     qnames:       String,
//     seq:          String, // of the top ref strand from start0 to end1
// }

// /// EncodingMetadata reports summary results of read encoding accumulated
// /// over all reads.
// pub struct EncodingMetadata {
//     pub n_reads:         usize,
//     pub n_unique_spans:  usize,
//     pub n_ref_bases:     usize,
//     pub n_read_bases:    usize,
//     pub n_match:         usize,
//     pub n_alt:           usize,
//     pub n_masked:        usize,
//     pub n_del:           usize,
//     pub n_ins:           usize,
//     pub n_reverse:       usize,
// }
// impl EncodingMetadata {
//     fn new(n_unique_spans: usize) -> Self {
//         EncodingMetadata {
//             n_reads:         0,
//             n_unique_spans:  n_unique_spans,
//             n_ref_bases:     0,
//             n_read_bases:    0,
//             n_match:         0,
//             n_alt:           0,
//             n_masked:        0,
//             n_del:           0,
//             n_ins:           0,
//             n_reverse:       0,
//         }
//     }
// }

// /// ReadEncodings aggregate one or more read encodings per ChromSpan for 
// /// subsequent sorting and printing.
// pub struct ReadEncodings(FxHashMap<ReFragment, Vec<ReadEncoding>>);
// impl ReadEncodings {
//     /// Create a new ReadEncodings object.
//     pub fn new() -> Self{
//         Self(FxHashMap::default())
//     }
//     /// Add a completed read encoding to the map.
//     pub fn insert(&mut self, chrom_index: u8, encoding: ReadEncoding) {
//         let re_fragment = ReFragment {
//             chrom_index: chrom_index,
//             start0:      encoding.start0,
//             end1:        encoding.end1,
//         };
//         self.0.entry(re_fragment)
//             .or_insert_with(Vec::new)
//             .push(encoding);          
//     }

//     /// Sort and write a vector of ReadEncodings to a temporary file.
//     pub fn write_sorted(
//         tool:   &SnvAnalysisTool,
//         worker: &mut SnvChromWorker,
//     ) -> EncodingMetadata {
//         let mut csv = OutputCsv::open_csv(
//             &worker.encodings_file_path, 
//             b'\t', 
//             false, 
//             Some(tool.n_cpu),
//         );
//         let mut fragments = worker.encodings.0.keys().filter_map(|s|{
//             let excluded  =  tool.exclusions.pos_in_region(&worker.chrom, s.start0 + 1);
//             let on_target = !tool.targets.has_data || 
//                                    tool.targets.pos_in_region(&worker.chrom, s.start0 + 1);
//             if !excluded && on_target { Some(s) } else { None }
//         }).collect::<Vec<_>>();
//         fragments.sort_unstable();
//         let mut md = EncodingMetadata::new(fragments.len());
//         for fragment in fragments {
//             let reads = &worker.encodings.0[&fragment];
//             let n_reads = reads.len();
//             let record = EncodingRecord {
//                 chrom_index: fragment.chrom_index,
//                 start0:      fragment.start0, // BED half-open coordinates
//                 end1:        fragment.end1,

//                 n_reads:     n_reads as u32,

//                 n_match:     reads.iter().map(|r| r.n_match).sum::<u32>(),
//                 n_alt:       reads.iter().map(|r| r.n_alt).sum::<u32>(),
//                 n_masked:    reads.iter().map(|r| r.n_masked).sum::<u32>(),
//                 n_del:       reads.iter().map(|r| r.n_del).sum::<u32>(),
//                 n_ins:       reads.iter().map(|r| r.n_ins).sum::<u32>(),
//                 n_reverse:   reads.iter().map(|r| r.strand as u32).sum::<u32>(),

//                 n_matches:   reads.iter().map(|r| r.n_match.to_string())
//                                 .collect::<Vec<String>>().join(","),
//                 n_alts:      reads.iter().map(|r| r.n_alt.to_string())
//                                 .collect::<Vec<String>>().join(","),
//                 n_maskeds:   reads.iter().map(|r| r.n_masked.to_string())
//                                 .collect::<Vec<String>>().join(","),
//                 n_dels:      reads.iter().map(|r| r.n_del.to_string())
//                                 .collect::<Vec<String>>().join(","),
//                 n_inss:      reads.iter().map(|r| r.n_ins.to_string())
//                                 .collect::<Vec<String>>().join(","),

//                 read_start0s: reads.iter().map(|r| r.read_start0.to_string())
//                                 .collect::<Vec<String>>().join(","),                   
//                 strands:      reads.iter().map(|r| r.strand.to_string())
//                                 .collect::<Vec<String>>().join(","),
//                 encodings:    reads.iter().map(|r| r.encoding.clone())
//                                 .collect::<Vec<String>>().join(","),
//                 insertions:   reads.iter().map(|r| r.insertion.join(""))
//                                 .collect::<Vec<String>>().join(","),

//                 sample_bits:  reads.iter().fold(0_u32, |acc, rd| acc | rd.sample_bit),
//                 sample_bitss: reads.iter().map(|r| r.sample_bit.to_string())
//                                 .collect::<Vec<String>>().join(","),
//                 qnames:       reads.iter().map(|r| r.qname.clone())
//                                 .collect::<Vec<String>>().join(","),

//                 seq: if n_reads >= 5 { // control file size
//                     tool.fa.view(
//                         worker.chrom_tid, 
//                         fragment.start0 as usize, 
//                         fragment.end1 as usize
//                     ).ok()
//                         .map(|v| v.to_string().to_ascii_uppercase())
//                         .unwrap_or_else(|| "NA".to_string())
//                 } else {
//                     "NA".to_string()
//                 }
//             };
//             csv.serialize(&record);
//             let fragment_len = (fragment.end1 - fragment.start0) as usize;
//             md.n_reads      += n_reads;
//             md.n_ref_bases  += fragment_len;
//             md.n_read_bases += n_reads * fragment_len;
//             md.n_match         += record.n_match as usize;
//             md.n_alt           += record.n_alt as usize;
//             md.n_masked        += record.n_masked as usize;
//             md.n_del           += record.n_del as usize;
//             md.n_ins           += record.n_ins as usize;
//             md.n_reverse       += record.n_reverse as usize;
//         }
//         md
//     }
// }


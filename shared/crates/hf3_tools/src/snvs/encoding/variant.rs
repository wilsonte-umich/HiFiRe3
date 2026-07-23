//! Support for creating an encoded variant-level file for visualization in the app.

// dependencies
// use std::str::from_utf8_unchecked;
// use rustc_hash::FxHashMap;
// use serde::Serialize;
// use mdi::OutputCsv;
// use super::{SnvChromWorker, SnvAnalysisTool};

// /// ReferenceVariantInstance holds the count of a specific ReferenceVariant 
// /// instance and the samples and reads that contributed to the count.
// struct ReferenceVariantInstance {
//     count:        usize,
//     sample_bits:  u32,
//     max_n_passes: u8,
//     n_allowed:    u8,
//     qnames:       Vec<String>,
// }
// impl ReferenceVariantInstance {
//     /// Create a new ReferenceVariantInstance instance with the specified count  
//     /// and sample bits.
//     pub fn new() -> Self {
//         ReferenceVariantInstance {
//             count:        0,
//             sample_bits:  0,
//             max_n_passes: 0,
//             n_allowed:    0,
//             qnames:       Vec::new(),
//         }
//     }
// }

// /// VariantMetadata reports summary results of variant calling and counting.
// pub struct VariantMetadata {
//     pub n_variants:       usize,
//     pub n_substitutions:  usize,
//     pub n_insertions:     usize,
//     pub n_deletions:      usize,
//     pub variant_count:    usize,
//     pub variant_coverage: usize,
//     pub is_simple_repeat: usize,
// }
// impl VariantMetadata {
//     fn new(n_variants: usize) -> Self {
//         VariantMetadata {
//             n_variants,
//             n_substitutions:  0,
//             n_insertions:     0,
//             n_deletions:      0,
//             variant_count:    0,
//             variant_coverage: 0,
//             is_simple_repeat: 0,
//         }
//     }
// }

// /// A VariantRecord is a specific SNV or indel, or a series of operations,
// /// observed at a single reference position on a chromosome, with a count of the 
// /// number of times the variant was observed at the position. Serializable for 
// /// writing to a tabix-compatible file.
// #[derive(Serialize)]
// struct VariantRecord<'a> {
//     chrom_index:  u8,
//     variant:      &'a ChromVariant, // specific variant observed at this position
//     count:        usize, // number of times this variant was observed at this position
//     coverage:     usize, // number of reads covering ref_pos0, including those without this variant
//     sample_bits:  u32,   // sample bits for this variant
//     n_samples:    u8,    // number of samples with this variant, derived from sample_bits
//     max_n_passes: u8,    // maximum (best) number of PacBio passes among reads with this variant
//     any_allowed:  u8,    // whether any read  with this variant passes duplex strand validation
//     all_allowed:  u8,    // whether all reads with this variant passes duplex strand validation
//     simple_repeat:u8,    // integer bool if this variant overlaps a simple sequence repeat
//     qnames:       String,// comma-delimited list of QNAMEs with this variant
// }

// /// ChromsVariants holds counts of specific ChromVariant instances 
// /// observed on a single chromosome.
// pub struct ChromVariants(FxHashMap<ChromVariant, ChromVariantObs>);
// impl ChromVariants {
//     /// Create a new ChromVariants instance with an empty hash map.
//     pub fn new() -> Self {
//         ChromVariants(FxHashMap::default())
//     }

//     /// Increment the count of a specific ChromVariant instance.
//     pub fn increment(
//         &mut self, 
//         ref_pos0:    u32, 
//         n_ref_bases: u32, 
//         alt_bases:   &str, 
//         sample_bit:  u32,
//         n_passes:    u8,
//         mut allowed: bool,
//         qname:       &[u8],
//     ) {
//         let variant = ChromVariant::new(ref_pos0, n_ref_bases, alt_bases);
//         let obs = self.0
//             .entry(variant.clone())
//             .or_insert_with(ChromVariantObs::new);
//         obs.count        += 1;
//         obs.sample_bits  |= sample_bit;
//         obs.max_n_passes =  obs.max_n_passes.max(n_passes);
//         obs.any_allowed  |= allowed;
//         obs.all_allowed  &= allowed;
//         unsafe {
//             obs.qnames.push(from_utf8_unchecked(qname).to_string());
//         }
//     }

//     /// Sort and write a vector of ChromVariant instances to a temporary file.
//     pub fn write_sorted(
//         tool:   &SnvAnalysisTool,
//         worker: &mut SnvChromWorker,
//     ) -> VariantMetadata {
//         let mut csv = OutputCsv::open_csv(
//             &worker.variants_file_path, 
//             b'\t', 
//             false, 
//             Some(tool.n_cpu),
//         );
//         let mut variants = worker.variants.0.keys().filter_map(|v|{
//             let excluded  =  tool.exclusions.pos_in_region(&worker.chrom, v.ref_pos0 + 1);
//             let on_target = !tool.targets.has_data || 
//                                    tool.targets.pos_in_region(&worker.chrom, v.ref_pos0 + 1);
//             if !excluded && on_target { Some(v.clone()) } else { None }
//         }).collect::<Vec<_>>();
//         variants.sort_unstable();
//         let mut md = VariantMetadata::new(variants.len());
//         for variant in variants {
//             let obs = &worker.variants.0[&variant];
//             let simple_repeat = worker.simple_repeats
//                 .binary_search(variant.ref_pos0, variant.n_ref_bases);
//             let record = VariantRecord {
//                 chrom_index:   worker.chrom_index,
//                 variant:       &variant,
//                 count:         obs.count,
//                 coverage:      worker.pileup.0[variant.ref_pos0 as usize].coverage(),
//                 sample_bits:   obs.sample_bits,
//                 n_samples:     obs.sample_bits.count_ones() as u8,
//                 max_n_passes:  obs.max_n_passes,
//                 any_allowed:   obs.any_allowed as u8,
//                 all_allowed:   obs.all_allowed as u8,
//                 simple_repeat: simple_repeat as u8,
//                 qnames:        obs.qnames.join(",")
//             };
//             csv.serialize(&record);
//             let alt_minus_ref = variant.alt_minus_ref();
//             if alt_minus_ref == 0 {
//                 md.n_substitutions += 1;
//             } else if alt_minus_ref > 0 {
//                 md.n_insertions += 1;
//             } else {
//                 md.n_deletions += 1;
//             }
//             md.variant_count    += record.count;
//             md.variant_coverage += record.coverage;
//             md.is_simple_repeat += simple_repeat as usize;
//         }
//         md
//     }
// }


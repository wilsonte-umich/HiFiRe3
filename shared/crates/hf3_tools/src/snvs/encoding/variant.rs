//! Support for creating encoded variant-level files for downstream use.

// dependencies
use rustc_hash::FxHashMap;
use serde::Serialize;
use mdi::OutputCsv;
use crate::snvs::*;

/// VariantZygosity lists the types of variant calls. Unlike encodings, a single
/// output file includes all variants calls.
#[derive(Clone, Copy, Serialize)]
#[repr(u8)]
pub enum VariantZygosity {
    Subclonal    = 0,
    Heterozygous = 1,
    Homozygous   = 2,
}

/// VariantInstances holds the read count of a specific Variant and the 
/// fragments, samples, and reads that contributed to the count.
pub struct VariantInstances {
    n_fragments:      usize, // the number of RE fragments that reported the variant (usually 1)
    n_matching_reads: usize,
    coverage:         usize,
    sample_bits:      SampleBits,
    max_avg_qual:     PhredQual,
    zygosity: VariantZygosity,
    vaf:      f64,
    qnames:   Vec<QName>,
}
impl VariantInstances {

    /// Create a new empty VariantInstances object.
    fn new() -> Self {
        VariantInstances {
            n_fragments:      0,
            n_matching_reads: 0,
            coverage:         0,
            sample_bits:      0,
            max_avg_qual:     0,
            zygosity: VariantZygosity::Subclonal,
            vaf:      0.0,
            qnames:   Vec::new(), // auto-allocate since many variants will have few reads
        }
    }
}

/// VariantMetadata reports summary results of variant calling and counting,
/// excluding variants in simple repeats.
pub struct VariantMetadata {
    pub n_variants:       usize,
    pub n_substitutions:  usize,
    pub n_insertions:     usize,
    pub n_deletions:      usize,
    pub n_homozygous:     usize,
    pub n_heterozygous:   usize,
    pub n_subclonal:      usize,
    pub variant_count:    usize,
    pub variant_coverage: usize,
}
impl VariantMetadata {
    fn new(n_variants: usize) -> Self {
        VariantMetadata {
            n_variants,
            n_substitutions:  0,
            n_insertions:     0,
            n_deletions:      0,
            n_homozygous:     0,
            n_heterozygous:   0,
            n_subclonal:      0,
            variant_count:    0,
            variant_coverage: 0,
        }
    }
}

/// A VariantRecord is a specific SNV or indel as written to file for 
/// downstream use. 
#[derive(Serialize)]
struct VariantRecord<'a> {
    chrom_index:      ChromIndex1,
    variant:          &'a Variant, // specific variant observed at this position
    n_fragments:      usize, // the number of RE fragments that reported the variant (usually 1)
    n_matching_reads: usize,
    coverage:         usize,
    sample_bits:      SampleBits,
    n_samples:        u32,
    zygosity:         VariantZygosity,
    vaf:              f64, // vaf only set on clonal
    max_avg_qual:     PhredQual, // max_avg_qual set on subclonal
    qnames: CommaDelimited, // comma-delimited list of QNAMEs with this variant
}

/// VariantsTally aggregates accumulated VariantInstances per Variant.
pub struct VariantsTally {
    pub tally: FxHashMap<Variant, VariantInstances>
}
impl VariantsTally {

    /// Create a new empty VariantsTally object. On variant tally is 
    /// instantiated per SnvChromWorker.
    pub fn new() -> Self {
        let mut tally = FxHashMap::default();
        tally.reserve(1_048_576);
        Self{ tally }
    }

    /// Update the instances tally of a specific clonal Variant derived from
    /// alignment to reference.
    pub fn add_clonal(
        &mut self, 
        variant: &Variant, 
        reads:   &[ReadInstance], // all ReFragment reads
        read_is: &[ReadIndex],    // indices into reads for reads with the variant
    ) {
        let instances = self.tally
            .entry(variant.clone())
            .or_insert_with(VariantInstances::new);
        instances.n_fragments      += 1;
        instances.n_matching_reads += read_is.len();
        instances.coverage         += reads.len();
        instances.vaf = instances.n_matching_reads as f64 / instances.coverage as f64;
        instances.zygosity = if instances.vaf > MAX_HETEROZYGOUS_ZYGOSITY {
            VariantZygosity::Homozygous
        } else {
            VariantZygosity::Heterozygous
        };
        for read_i in read_is {
            let read = &reads[*read_i];
            instances.sample_bits |= read.sample_bit;
            instances.qnames.push(read.qname.clone());
        }
    }

    /// Update the instances tally of a specific subclonal Variant derived from
    /// alignment to a haplotype consensus.
    pub fn add_subclonal(
        &mut self, 
        variant: &Variant, 
        reads:   &[ReadInstance], // all ReFragment reads
        read_is: &[ReadIndex],    // read indices for the haplotype being processed
        read_js: &[ReadIndex],    // indices into read_is for reads with the variant
        max_avg_qual: PhredQual,
    ) {
        let instances = self.tally
            .entry(variant.clone())
            .or_insert_with(VariantInstances::new);
        instances.n_fragments      += 1;
        instances.n_matching_reads += read_js.len();
        instances.coverage         += read_is.len();
        instances.max_avg_qual = instances.max_avg_qual.max(max_avg_qual);
        for read_j in read_js {
            let read = &reads[read_is[*read_j]];
            instances.sample_bits |= read.sample_bit;
            instances.qnames.push(read.qname.clone());
        }
    }

    /// Sort and write a set of VariantInstances to a temporary file for the
    /// working chromosome.
    pub fn write_sorted(
        tool:   &SnvAnalysisTool,
        worker: &mut SnvChromWorker,
    ) -> VariantMetadata {
        let mut csv = OutputCsv::open_csv(
            &worker.variants_file_path, 
            b'\t', 
            false, 
            Some(tool.n_cpu),
        );
        let mut variants = worker.variant_tally.tally.keys()
            .filter_map(|v|{
                let excluded  =  tool.exclusions.pos_in_region(&worker.chrom, v.tgt_pos0 + 1);
                let on_target = !tool.targets.has_data || 
                                       tool.targets.pos_in_region(&worker.chrom, v.tgt_pos0 + 1);
                if !excluded && on_target { Some(v.clone()) } else { None }
            }).collect::<Vec<_>>();
        variants.sort_unstable();
        let mut md = VariantMetadata::new(variants.len());
        for variant in variants {
            let instances = &worker.variant_tally.tally[&variant];
            let record = VariantRecord {
                chrom_index:   worker.chrom_index,
                variant:       &variant,
                n_fragments:      instances.n_fragments,
                n_matching_reads: instances.n_matching_reads,
                coverage:      instances.coverage,
                sample_bits:   instances.sample_bits,
                n_samples:     instances.sample_bits.count_ones(),
                zygosity:      instances.zygosity,
                vaf:           instances.vaf,
                max_avg_qual:  instances.max_avg_qual,
                // any_allowed:   instances.any_allowed as u8,
                // all_allowed:   instances.all_allowed as u8,
                qnames:   instances.qnames.join(",")
            };
            csv.serialize(&record);

            let alt_minus_ref = variant.alt_minus_ref();
            if alt_minus_ref == 0 {
                md.n_substitutions += 1;
            } else if alt_minus_ref > 0 {
                md.n_insertions += 1;
            } else {
                md.n_deletions += 1;
            }
            match record.zygosity {
                VariantZygosity::Subclonal => {
                    md.n_subclonal += 1
                },
                VariantZygosity::Heterozygous => {
                    md.n_heterozygous += 1
                },
                VariantZygosity::Homozygous => {
                    md.n_homozygous += 1
                }
            }
            md.variant_count    += record.n_matching_reads;
            md.variant_coverage += record.coverage; 
        }
        md
    }
}

//! Support for building haplotype consensuses from RE fragments and using them
//! to call rare variants without reference alignment bias.

// dependencies
use std::iter::repeat_n;
use rustc_hash::FxHashMap;
// use serde::Serialize;
use genomex::sequence::AlignmentStatus;
use super::*;

// constants
const MIN_HAPLOTYPE_READS:        usize = 3;  // a haplotype is defined as having >=3 matching reads
const MIN_HETEROZYGOUS_INSTANCES: usize = 8;  // allow fewer reads if heterozygosity is detected
const MIN_HOMOZYGOUS_INSTANCES:   usize = 15; // require more homozygous reads to minimize false negative heterozygosity

impl SnvChromWorker {

    /// Build homozygous or heterozygous consensuses for RE fragment haplotypes.
    pub fn build_consensuses(
        &mut self,
        fragment_reads: FragmentReads,
    ){
        for re_fragment in fragment_reads.instances.keys() {

            // require a minimum fragment coverage to continue
            let reads = fragment_reads.instances.get(&re_fragment).unwrap();
            let n_reads = reads.len();
            if n_reads < MIN_HETEROZYGOUS_INSTANCES { continue; }

            // collect all variants over all reads into a variant by read map
            let mut ref_vars = FragmentReferenceVariants::new();
            reads.iter().enumerate().for_each(|(read_i, read)|{
                self.fill_fragment_variants(&mut ref_vars, read_i, read.start0, &read.cs);
            });

            // restrict attention to only heterozygous variants
            // disregard homozygous variants and rare variants for now
            let rvars: Vec<_> = ref_vars.variant_map.keys().cloned().collect();
            for rvar in rvars{
                let vmap =  &ref_vars.variant_map[&rvar];
                let vaf = vmap.n_matching_reads as f64 / n_reads as f64;
                let zyg = ((vaf - 0.5).abs() - 0.5).abs();
                if vmap.n_matching_reads < MIN_HETEROZYGOUS_INSTANCES || zyg < 0.2 {
                    ref_vars.variant_map.remove(&rvar);
                }
            }

            // require more reads to continue with a homozygous RE fragment
            let n_heterozygous_variants = ref_vars.variant_map.len();
            if n_heterozygous_variants == 0 {
                if n_reads < MIN_HOMOZYGOUS_INSTANCES { continue; }
                let hap1_read_is: Vec<ReadIndex> = (0..n_reads).map(|i| i). collect();
                let hap1_seq = self.build_haplotype_consensus(&reads, &hap1_read_is);
                // fill_haplo
                continue;
            }

            // normalize non-reference variants to consistent haplotype numbers
            // given that variant SNPs may be on opposing haplotypes
            // reads that cannot be unamiguosly assigned are haplotype 0 at this point
            let het_vars: Vec<_> = ref_vars.variant_map.keys().cloned().collect();
            let index_matches = ref_vars.variant_map[&het_vars[0]].read_map.clone();
            let mut read_haplotypes: Vec<u8> = if n_heterozygous_variants > 1 { 
                for het_var in het_vars[1..n_heterozygous_variants].iter(){
                    let vmap =  ref_vars.variant_map.get_mut(&het_var).unwrap();
                    let match_score = vmap.read_map.iter().zip(&index_matches).map(|x|{
                        if *x.0 == *x.1 { 1 } else { -1 }
                    }).sum::<i32>();
                    if match_score < 0 {
                        for b in vmap.read_map.iter_mut() { *b = !*b; }
                    }
                }
                (0..n_reads).map(|read_i|{
                    let rmap: Vec<bool> = het_vars.iter().map(|het_var|{
                        let vmap =  ref_vars.variant_map.get(&het_var).unwrap();
                        vmap.read_map[read_i]
                    }).collect();
                    if rmap.windows(2).all(|w| w[0] == w[1]) {
                        rmap[0] as u8 + 1
                    } else {
                        0 // a read wasn't always assigned to the same haplotype
                    }
                }).collect()
            } else {
                index_matches.into_iter().map(|b| b as u8 + 1).collect()
            };

            // abort if we failed to resolve a fragment into haplotype reads for consensus
            let mut hap1_read_is: Vec<ReadIndex> = read_haplotypes.iter().enumerate().filter_map(|(read_i, hap)|{
                if *hap == 1 { Some(read_i) } else { None }
            }).collect();
            let mut hap2_read_is: Vec<ReadIndex> = read_haplotypes.iter().enumerate().filter_map(|(read_i, hap)|{
                if *hap == 2 { Some(read_i) } else { None }
            }).collect();
            if hap1_read_is.len() < MIN_HAPLOTYPE_READS || 
               hap2_read_is.len() < MIN_HAPLOTYPE_READS { 
                continue;
            }

            // build consensus sequences of each haplotype
            let hap1_seq = self.build_haplotype_consensus(&reads, &hap1_read_is);
            let hap2_seq = self.build_haplotype_consensus(&reads, &hap2_read_is);

            // assign ambiguous reads to a haplotype
            let ambiguous_read_is: Vec<ReadIndex> = read_haplotypes.iter().enumerate().filter_map(|(read_i, hap)|{
                if *hap == 0 { Some(read_i) } else { None }
            }).collect();
            if ambiguous_read_is.len() > 0 {
                ambiguous_read_is.iter().for_each(|read_i| {
                    read_haplotypes[*read_i] = self.assign_ambiguous_to_haplotype(&reads, *read_i, &hap1_seq, &hap2_seq);
                });    
                hap1_read_is = read_haplotypes.iter().enumerate().filter_map(|(read_i, hap)|{
                    if *hap == 1 { Some(read_i) } else { None }
                }).collect();
                hap2_read_is = read_haplotypes.iter().enumerate().filter_map(|(read_i, hap)|{
                    if *hap == 2 { Some(read_i) } else { None }
                }).collect();            
            }

            // align each read to its haplotype consensus to call final rare variants
            let mut hap_vars = FragmentHaplotypeVariants::new();
            // self.fill_haplotype_variants(&read,);
            


        }
    }

    /// Use the unambiguous haplotype reads to build a haplotype consensus sequence.
    fn build_haplotype_consensus(
        &mut self,
        reads:   &[ReadInstance],
        read_is: &[ReadIndex],
    ) -> String {
        let seq0 = reads[read_is[0]].get_top_strand_seq();
        let n_bases = seq0.len();
        let n_reads = reads.len();
        let mut map: Vec<Vec<String>> = Vec::with_capacity(n_reads);
        map.push(seq0.split("").map(|c| c.to_string()).collect());
        for read_j in 1..n_reads{
            let aln = self.aligner.align(
                &reads[read_is[read_j]].get_top_strand_seq(), 
                &seq0, 
                None, 
                false,
            );
            if aln.status == AlignmentStatus::AlignmentFound {
                let mut qry_on_tgt = Vec::with_capacity(n_bases);
                if aln.tgt_start0 > 0 {
                    qry_on_tgt.extend(repeat_n("N".to_string(), aln.tgt_start0));
                }
                qry_on_tgt.extend(aln.qry_on_tgt);
                if aln.tgt_end0 < n_bases - 1 {
                    qry_on_tgt.extend(repeat_n("N".to_string(), n_bases - 1 - aln.tgt_end0));
                }
                map.push(qry_on_tgt);                
            }
        }
        let n_reads = map.len();
        (0..n_bases)
            .map(|base_i|{
                (0..n_reads)
                    .map(|read_j| &map[read_j][base_i])
                    .fold(FxHashMap::default(), |mut counts, base| {
                        *counts.entry(base).or_insert(0) += 1;
                        counts
                    })
                    .into_iter()
                    .max_by(|a, b| a.1.cmp(&b.1))
                    .map(|(base, _)| base)
                    .unwrap()
            })
            .filter_map(|base|{
                if base == "-" { None} else { Some(base.to_owned()) }
            })
            .collect::<Vec<String>>()
            .join("")
    }

    /// Use the unambiguous haplotype reads to build a haplotype consensus sequence.
    fn assign_ambiguous_to_haplotype(
        &mut self,
        reads:  &[ReadInstance],
        read_i: ReadIndex,
        hap1_seq: &str,
        hap2_seq: &str,
    ) -> u8 {
        let seq = &reads[read_i].get_top_strand_seq();
        let aln1 = self.aligner.align(
            &seq, 
            hap1_seq, 
            None, 
            false,
        );
        let aln2 = self.aligner.align(
            &seq, 
            hap2_seq, 
            None, 
            false,
        );
        if aln1.score > aln2.score { 1 } else { 2 }
    }

    fn fill_haplotype_variants(
        &self,
        fvars: &mut FragmentReferenceVariants,
        read_i: ReadIndex,
        aln_start0: u32,
        cs: &str,
    ) { 

    }

    /// Align one read sequence to the haplotype consensus.
    fn align_to_haplotype_consensus(
        &mut self,
        reads:   &[ReadInstance],
        read_is: &[ReadIndex],
        hap_seq: &str,
    ) {
        for read_i in read_is {
            let aln = self.aligner.align(
                &reads[*read_i].get_top_strand_seq(), 
                hap_seq, 
                None, 
                false,
            );

        }
    }

}



        //     // build a map of all variant read bases on the reference sequence from cs tags
        //     let mut reference_positions = ReferencePositions(
        //         vec![ReferencePosition::new(); n_ref_bases]
        //     );
        //     reads.iter().enumerate().for_each(|(read_i, read)|{
        //         self.set_read_on_ref(read.start0, re_fragment.start0, &read.cs);
        //         for ref_pos_i in 0..n_ref_bases{
        //             let Some(Some(bases)) = self.read_on_ref.get(ref_pos_i) else{ continue; };
        //             let ref_pos = &mut reference_positions.0[ref_pos_i];
        //             ref_pos.n_allowed_reads += 1;
        //             if bases != "=" { 
        //                 ref_pos.insert(bases, read_i); 
        //             }
        //         }
        //     });

        //     // let mut heterozygous_matrix: Vec<Vec<bool>> = Vec::new();
        //     // let mut heterozygous_bases: Vec<String> = Vec::new();
        //     // let mut sum_zygosity = 0.0;
        //     // reference_positions.0.iter().for_each(|ref_pos|{
        //     //     for (variant_bases, read_is) in &ref_pos.allowed_variants {
        //     //         let n_instances = read_is.len();
        //     //         let vaf = n_instances as f64 / ref_pos.n_allowed_reads as f64;
        //     //         let zyg = ((vaf - 0.5).abs() - 0.5).abs();
        //     //         if n_instances >= 3 && zyg >= 0.2 {
        //     //             let mut x = vec![false; n_reads];
        //     //             read_is.iter().for_each(|read_i| x[*read_i] = true );
        //     //             heterozygous_matrix.push(x);
        //     //             heterozygous_bases.push(variant_bases.clone());
        //     //             sum_zygosity += zyg;
        //     //         }
        //     //     }
        //     // });
            
        //     // let n_heterozygous = heterozygous_matrix.len();
        //     // if n_heterozygous == 0 {
        //     //     // process homozygous
        //     // } else if n_heterozygous == 1 {
        //     //     // use haplotypes as set by one variant
        //     // } else {
        //     //     let mean_zygosity = sum_zygosity / n_heterozygous as f64;
        //     //     // use votes to find one consistent reference read for each haplotype
        //     //     // not all reads need to be consistent in each haplotype, just one per haplotype

        //     // }

        //     // // finish constructing the many-to-many relationship between reads and variants
        //     // let mut variant_instances: FxHashMap::<ReferenceVariant, Vec<usize>> = FxHashMap::default();
        //     // let mut variant_allowed_instances: FxHashMap::<ReferenceVariant, Vec<usize>> = FxHashMap::default();
        //     // reads.iter().enumerate().for_each(|(read_i, read)|{
        //     //     read.reference_variants.iter().for_each(|(var, allowed)|{
        //     //         variant_instances.entry(var.clone())
        //     //             .or_insert_with(Vec::new)
        //     //             .push(read_i);  
        //     //         if *allowed {
        //     //             variant_allowed_instances.entry(var.clone())
        //     //                 .or_insert_with(Vec::new)
        //     //                 .push(read_i); 
        //     //         } 
        //     //     });
        //     // });

        //     // // handle rare cases where a fragment called no variants in any read
        //     // // must be homozygous reference
        //     // if variant_instances.len() == 0 {
        //     //     ReFragment::build_consensus(reads);
        //     //     continue;
        //     // }

        //     // let mut allowed_variant_counts: FxHashSet<usize> = FxHashSet::default();
        //     // let mut max_allowed_variant_count = 0;
        //     // variant_allowed_instances.values().for_each(|read_is| {
        //     //     let n_allowed_instances = read_is.len();
        //     //     allowed_variant_counts.insert(n_allowed_instances);
        //     //     max_allowed_variant_count = max_allowed_variant_count.max(n_allowed_instances);
        //     // });
        //     // let mut max_other_count = 0;
        //     // variant_instances.values().for_each(|read_is| {
        //     //     let n_other_instances = n_reads - read_is.len(); // usually reference values
        //     //     max_other_count = max_other_count.max(n_other_instances);
        //     // });


        //     // if max_allowed_variant_count < MIN_HAPLOTYPE_READS || 
        //     //    max_other_count < MIN_HAPLOTYPE_READS {
        //     //     ReFragment::build_consensus(reads);
        //     //     continue;
        //     // }

        //     // // use heterozygous variants to cluster reads into haplotypes
        //     // let variants = fragment_variants.0.get(re_fragment);
        //     // if let Some(variants) = variants {

        //     //     // variants.iter().map(|var|{
        //     //     //     var.qnames
        //     //     // })
        //     //     let mut variant_instances = FxHashMap::<ReferenceVariant, Vec<usize>>::default();
        //     //     // for read in reads {
        //     //     //     for variant in read.reference_variants {

        //     //     //     }
        //     //     // }

        //     //     reads.iter().enumerate().for_each(|(read_i, read)|{
        //     //         read.reference_variants.iter().for_each(|var|{
        //     //             variant_instances.entry(var.clone())
        //     //                 .or_insert_with(Vec::new)
        //     //                 .push(read_i);  
        //     //         });
        //     //     });


                

        //     // // handle rare cases where a fragment called no variants in any read
        //     // // must be homozygous reference
        //     // } else {
        //     //     ReFragment::build_consensus(reads);
        //     // }


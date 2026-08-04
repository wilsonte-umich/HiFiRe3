//! Support for aligning reads to the haplotype consensus for analyze_reads.

// imports
use std::cmp::Ordering;
use minimap2::{Aligner as Minimap2, Built};
use crate::snvs::*;
use super::*;

// constants
const ONE_THIRD:  f64 = 1.0 / 3.0;
const TWO_THIRDS: f64 = 2.0 / 3.0;

pub struct ReadAssignment {
    pub haplotype:    Haplotype,
    pub cs:           String,
    pub query_start:  i32,
    pub target_start: i32,
} 

impl SnvChromWorker {

    /// Assign one read to each haplotype consensus to resolve its haplotype.
    pub(super) fn assign_read_to_haplotype(
        &mut self,
        read:     &ReadInstance,
        read_i:   ReadIndex,
        mm2_hap1: &Minimap2<Built>,
        mm2_hap2: &Minimap2<Built>,
    ) -> Option<ReadAssignment> {
        // if self.show_debug {
        //     eprintln!("");
        //     eprintln!("read_i {read_i}");
        // }
        // collect haplotype votes for this read by comparing 
        // read_on_ref tracking variants to ref_on_hap variants
        self.reset_hap_votes();
        for variant in &self.tracking_variants {
            let vmap = self.frag_vars.variant_map.get_mut(&variant).unwrap();
            let haplotype_matches = (
                self.hap_vars[&Haplotype::Haplotype1].contains(variant),
                self.hap_vars[&Haplotype::Haplotype2].contains(variant)
            );
            // if self.show_debug &&
            //    (vmap.read_map[read_i].has_var() || 
            //     vmap.read_map[read_i].is_informative) {
            //     eprintln!("variant {:?}", variant);
            //     eprintln!("haplotype_matches {:?}", haplotype_matches);
            // }
            if vmap.read_map[read_i].has_var() {
                match haplotype_matches {
                    m if m == (true, false) => {
                        *self.hap_votes.get_mut(&Haplotype::Haplotype1).unwrap() += 1;
                    },
                    m if m == (false, true) => {
                        *self.hap_votes.get_mut(&Haplotype::Haplotype2).unwrap() += 1;
                    },
                    // homozygous (true, true) and subclonal (false, false) don't vote
                    _ => {} 
                }
            } else if vmap.read_map[read_i].is_informative {
                match haplotype_matches {
                    m if m == (false, true) => {
                        *self.hap_votes.get_mut(&Haplotype::Haplotype1).unwrap() += 1;
                    },
                    m if m == (true, false) => {
                        *self.hap_votes.get_mut(&Haplotype::Haplotype2).unwrap() += 1;
                    },
                    _ => {}
                }
            }
            // non-informative variant no-calls don't vote
        }   
        // if self.show_debug {
        //     eprintln!("self.hap_votes {:?}", self.hap_votes);
        // }

        // read assignments are resolved by majority vote across informative 
        // variants when the bias is sufficient, otherwise they are discarded
        let n_hap1_votes = self.hap_votes[&Haplotype::Haplotype1] as f64;
        let n_hap2_votes = self.hap_votes[&Haplotype::Haplotype2] as f64;
        let frac_hap1 = n_hap1_votes / (n_hap1_votes + n_hap2_votes);
        if frac_hap1 >= TWO_THIRDS  { // 1 of 1, 2 of 3, 3 of 4, 4 of 5, 4 of 6, 5 of 7, 6 of 8, etc.
            Self::get_read_haplotype(read, mm2_hap1, Haplotype::Haplotype1)
        } else if frac_hap1 <= ONE_THIRD {
            Self::get_read_haplotype(read, mm2_hap2, Haplotype::Haplotype2)
        } else {
            // reads with no votes end here since frac_hap1 is NaN, which compares false
            // also ambiguous reads with 1 of 2, 2 of 4, 3 of 5, 5 of 8, etc.
            None 
        }
    }

    /// Get the read assignment for an unambiguously assigned read haplotype.
    fn get_read_haplotype(
        read:      &ReadInstance,
        mm2_hap:   &Minimap2<Built>,
        haplotype: Haplotype,
    ) -> Option<ReadAssignment> {
        let read_on_hap = &mm2_hap.map(
            &read.seq_bytes, 
            true, // cs not used, but need full score calculation
            false, 
            None, 
            Some(&MM_F_NO_PRINT_2ND), 
            None
        ).expect("Minimap2 failed at read_on_hap")[0];
        let Some(aln) = &read_on_hap.alignment else { return None; };
        let Some(cs) = &aln.cs else { return None; }; // don't expect failures
        Some(ReadAssignment{
            haplotype, 
            cs: cs.clone(), 
            query_start:  read_on_hap.query_start,
            target_start: read_on_hap.target_start,
        })
    }

    /// Get the read assignment for read with conflicting haplotype votes
    /// by alignment score. Deprecated: too prone to wrong conclusions.
    fn _get_read_haplotype_by_score(
        &mut self,
        read:     &ReadInstance,
        mm2_hap1: &Minimap2<Built>,
        mm2_hap2: &Minimap2<Built>,
    ) -> Option<ReadAssignment> {
        let read_on_hap1 = &mm2_hap1.map(
            &read.seq_bytes, 
            true, // cs not used, but need full score calculation
            false, 
            None, 
            Some(&MM_F_NO_PRINT_2ND), 
            None
        ).expect("Minimap2 failed at read_on_hap1")[0];
        let read_on_hap2 = &mm2_hap2.map(
            &read.seq_bytes, 
            true, // cs not used, but need full score calculation
            false, 
            None, 
            Some(&MM_F_NO_PRINT_2ND), 
            None
        ).expect("Minimap2 failed at read_on_hap2")[0];
        let Some(aln1) = &read_on_hap1.alignment else { return None; };
        let Some(aln2) = &read_on_hap2.alignment else { return None; };
        let Some(score1) = &aln1.alignment_score else { return None; };
        let Some(score2) = &aln2.alignment_score else { return None; };
        let Some(cs1) = &aln1.cs else { return None; }; // don't expect failures
        let Some(cs2) = &aln2.cs else { return None; };
        match score1.cmp(score2) {
            Ordering::Greater => Some(ReadAssignment{
                haplotype: Haplotype::Haplotype1, 
                cs: cs1.clone(), 
                query_start:  read_on_hap1.query_start,
                target_start: read_on_hap1.target_start,
            }),
            Ordering::Less => Some(ReadAssignment{
                haplotype: Haplotype::Haplotype2, 
                cs: cs2.clone(), 
                query_start:  read_on_hap2.query_start,
                target_start: read_on_hap2.target_start,
            }),
            // rare reads where we truly can't determine the haplotype are dropped
            Ordering::Equal => None,
        }

    }

    /// Align all haplotype read sequences to their consensus to call subclonal
    /// variants relative to that consensus.
    pub(super) fn align_to_haplotype_consensus(
        &mut self,
        reads_on_haplotype: &mut FragmentHaplotypes,
        re_fragment:  &ReFragment,
        haplotype:    Haplotype,
        ref_pos0_map: &mut Vec<ChromPos0>,
        reads:        &[ReadInstance],
        read_is:      Vec<ReadIndex>,
        read_masks:   &Vec<Vec<DdMaskType>>,
        mm2_hap:      Option<Minimap2<Built>>,
        assignments:  Option<Vec<ReadAssignment>>,
    ) {

        // align each haplotype read to its consensus if not already done
        let n_reads = reads.len();
        let n_haplotype_reads = read_is.len();
        self.frag_vars.reset(n_haplotype_reads);
        for read_j in 0..n_haplotype_reads {
            let read_i = read_is[read_j];
            let read = &reads[read_i];
            if let Some(assignments) = &assignments {
                let assignment = &assignments[read_j];
                self.encoding.prepare_read_on_hap(
                    re_fragment, read, assignment.target_start as usize
                );
                self.process_cs_tag(
                    reads_on_haplotype, re_fragment, haplotype, 
                    Some((assignment.query_start, assignment.target_start)), 
                    Some(&assignment.cs), ref_pos0_map,
                    // not read_i; read_j indexes into frag_vars
                    read_j, read, &read_masks[read_i] 
                );  
            } else { // this path is used by fully homozygous fragments
                let mm2_hap = mm2_hap.as_ref().unwrap();
                let read_on_hap = &mm2_hap.map(
                    &read.seq_bytes, 
                    true, 
                    false, 
                    None, 
                    Some(&MM_F_NO_PRINT_2ND), 
                    None
                ).expect("Minimap2 failed at read_on_hap")[0];
                let Some(aln) = &read_on_hap.alignment else { return; };
                let Some(cs) = &aln.cs else { return; }; // never expected to fail
                self.encoding.prepare_read_on_hap(
                    re_fragment, read, read_on_hap.target_start as usize
                );
                self.process_cs_tag(
                    reads_on_haplotype, re_fragment, haplotype, 
                    Some((read_on_hap.query_start, read_on_hap.target_start)), 
                    Some(cs), ref_pos0_map,
                    // not read_i; read_j indexes into frag_vars
                    read_j, read, &read_masks[read_i] 
                );     
            }
            reads_on_haplotype.insert_encoding(
                re_fragment, haplotype, self.encoding.clone()
            ); 
        }

        // call subclonal variants aggregated over all haplotype reads (might be none)
        for variant in self.frag_vars.variant_map.keys(){
            let vmap =  &self.frag_vars.variant_map[&variant];
            let read_js: Vec<ReadIndex> = vmap.read_map.iter()
                .enumerate()
                .filter_map(|(read_j, r)|{
                    if r.has_var() { Some(read_j) } else { None }
                }).collect();
            let max_avg_qual = read_js.iter()
                .map(|read_j| {
                    let read_i = read_is[*read_j];
                    let avg_qual= vmap.read_map[*read_j].avg_qual;
                    self.variant_reads_tally.add_subclonal_variant(
                        &reads[read_i], re_fragment, &haplotype, 
                        variant, avg_qual, 
                        n_haplotype_reads, n_reads
                    );
                    avg_qual
                })
                .max()
                .unwrap_or_default();
            self.variant_tally.add_subclonal(
                &variant, reads, &read_is, &read_js, 
                max_avg_qual
            );
        }
    }
}

//! Support for building haplotype consensuses from RE fragments and using them
//! to call rare variants without reference alignment bias.
//! 
//! Clonal variants are called from recurring variants seen in prior minimap2
//! alignments against reference. Subclonal variants are called from new
//! alignments against haplotype consensus(es) generated here.
//! 
//! `analyze_reads` in this script consumes the most time in `analyze SNVs` 
//! action execution.

// imports
use std::iter::repeat_n;
use rustc_hash::FxHashMap;
use super::*;

// constants
const MAX_EXPECTED_READ_LEN:  usize = 10000;
const MIN_HAPLOTYPE_READS:    usize = 3;  // a haplotype must have >=3 matching reads
const MIN_HETEROZYGOUS_READS: usize = 8;  // allow fewer total reads when heterozygous
const MIN_HOMOZYGOUS_READS:   usize = 15; // require more homozygous reads to minimize false homozygosity
const MM_F_NO_PRINT_2ND: [u64; 1] = [16384]; // minimap2 flag to suppress return of secondary alignments

impl SnvChromWorker {

    /// Build homozygous or heterozygous consensuses for RE fragment haplotypes,
    /// and use them to call clonal and subclonal variants.
    pub fn analyze_reads(
        &mut self,
        tool: &SnvAnalysisTool,
        fragment_reads: FragmentReads,
        haplotype_consensuses: &mut HaplotypeConsensuses,
    ){
        // allocate recycled objects used in fragment read parsing
        let mut frag_vars = FragmentVariants::new();
        let mut encoding = AlignmentEncoding::new();
        let mut seq0_bases: Vec<String> = Vec::with_capacity(MAX_EXPECTED_READ_LEN);
        let mut cs_map: Vec<FxHashMap<String, u8>> = Vec::with_capacity(MAX_EXPECTED_READ_LEN);
        let mut op_val: String = String::with_capacity(128);
        let mut ref_pos0_map: Vec<ChromPos0> = Vec::with_capacity(MAX_EXPECTED_READ_LEN);
        let mut hap_vs_ref = String::with_capacity(128);

        // process all observed ReFragents one at a time
        for re_fragment in fragment_reads.instances.keys() {

            // require a minimum RE fragment coverage to continue
            let reads = fragment_reads.instances.get(&re_fragment).unwrap();
            let n_reads = reads.len();
            if n_reads < MIN_HETEROZYGOUS_READS { continue; }

            // cache the RE fragment reference sequence for use below and printing
            haplotype_consensuses.insert_reference(tool, self, re_fragment);

            // collect all variants over all reads into a variant-by-read matrix
            frag_vars.reset(n_reads);
            reads.iter().enumerate().for_each(|(read_i, read)|{
                encoding.prepare_read_on_ref(re_fragment, read);
                self.process_cs_tag(
                    Haplotype::Unspecified, None, None, None,
                    &mut frag_vars, &mut encoding, 
                    read_i, read
                );
                self.reads_on_reference.insert(
                    re_fragment, 
                    Haplotype::Unspecified,
                    encoding.clone()
                );
            });

            // commit homozyogus and heterozygous variants to outputs
            // then restrict attention to only heterozygous variants for now
            let ref_vars: Vec<_> = frag_vars.variant_map.keys().cloned().collect();
            for ref_var in ref_vars{
                let vmap =  &frag_vars.variant_map[&ref_var];
                if vmap.n_matching_reads >= MIN_HAPLOTYPE_READS { // variant is at least heterozygous
                    let read_is: Vec<ReadIndex> = vmap.read_map.iter()
                        .enumerate()
                        .filter_map(|(read_i, r)|{
                            if r.has_var() { Some(read_i) } else { None }
                        }).collect();
                    self.variant_tally.add_clonal(
                        &ref_var, reads, &read_is
                    );
                    let vaf = vmap.n_matching_reads as f64 / n_reads as f64;
                    if vaf > MAX_HETEROZYGOUS_ZYGOSITY {
                        frag_vars.variant_map.remove(&ref_var); // done with homozygous
                    }
                } else {
                    frag_vars.variant_map.remove(&ref_var); // disregard subclonal for now
                }
            }

            // short-circuit to final processing homozygous RE fragments
            let n_heterozygous_variants = frag_vars.variant_map.len();
            if n_heterozygous_variants == 0 {
                // require more reads to continue with a homozygous RE fragment
                if n_reads < MIN_HOMOZYGOUS_READS { continue; }
                let hap3_read_is: Vec<ReadIndex> = (0..n_reads).map(|i| i). collect();
                let hap3_seq = self.build_haplotype_consensus(
                    &mut seq0_bases, &mut cs_map, &mut op_val,
                    &reads, &hap3_read_is
                );
                let Some(hap3_seq) = hap3_seq else { continue; };
                self.align_to_haplotype_consensus(
                    &mut op_val, &mut ref_pos0_map, &mut hap_vs_ref, 
                    haplotype_consensuses, &mut frag_vars, &mut encoding, 
                    re_fragment, reads, 
                    hap3_read_is, hap3_seq, Haplotype::Homozygous
                );
                continue;
            }

            // normalize non-reference variants to consistent haplotype numbers
            // given that different variants may be on opposing haplotypes
            // reads that cannot be unambiguously assigned are haplotype 0 at this point
            // this process might change mutable `hap0` but does not change immutable `has_var`
            let het_vars: Vec<_> = frag_vars.variant_map.keys().cloned().collect();
            let index_matches = frag_vars.variant_map[&het_vars[0]].read_map.clone();
            let mut read_haplotypes: Vec<u8> = if n_heterozygous_variants > 1 { 
                // re-orient variant read map hap0s to index_matches, i.e., to the first variant
                for het_var in het_vars[1..n_heterozygous_variants].iter(){
                    let vmap =  frag_vars.variant_map.get_mut(&het_var).unwrap();
                    let match_score = vmap.read_map.iter()
                        .zip(&index_matches)
                        .map(|(r1, r2)|{
                            if r1.has_var() == r2.has_var() { 1 } else { -1 }
                        }).sum::<i32>();
                    if match_score < 0 {
                        for r in vmap.read_map.iter_mut() { 
                            r.hap0 = !r.hap0; 
                        }
                    }
                }
                // assign each unambiguous read to haplotype 1 or 2; 0 if ambiguous
                (0..n_reads).map(|read_i|{
                    let rmap: Vec<ReadMapEntry> = het_vars.iter()
                        .map(|v|{
                            let vmap =  frag_vars.variant_map.get(&v).unwrap();
                            vmap.read_map[read_i]
                        }).collect();
                    let all_haps_match = rmap
                        .windows(2)
                        .all(|rs| rs[0].hap0 == rs[1].hap0);
                    if all_haps_match {
                        rmap[0].hap0 as u8 + 1
                    } else {
                        0 // a read wasn't always assigned to the same haplotype (for now)
                    }
                }).collect()
            } else {
                // assign single heterozygous variants as haplotype2 
                index_matches.into_iter()
                    .map(|r| r.hap0 as u8 + 1).collect()
            };

            // abort if we failed to resolve a fragment into haplotype reads for consensus
            let mut hap1_read_is: Vec<ReadIndex> = read_haplotypes.iter().enumerate()
                .filter_map(|(read_i, hap)|{
                    if *hap == 1 { Some(read_i) } else { None }
                }).collect();
            let mut hap2_read_is: Vec<ReadIndex> = read_haplotypes.iter().enumerate()
                .filter_map(|(read_i, hap)|{
                    if *hap == 2 { Some(read_i) } else { None }
                }).collect();
            if hap1_read_is.len() < MIN_HAPLOTYPE_READS || 
               hap2_read_is.len() < MIN_HAPLOTYPE_READS { 
                continue;
            }

            // build consensus sequences of each haplotype using unambiguous reads
            let hap1_seq = self.build_haplotype_consensus(
                &mut seq0_bases, &mut cs_map, &mut op_val,
                &reads, &hap1_read_is
            );
            let hap2_seq = self.build_haplotype_consensus(
                &mut seq0_bases, &mut cs_map, &mut op_val,
                &reads, &hap2_read_is
            );
            let Some(hap1_seq) = hap1_seq else { continue; };
            let Some(hap2_seq) = hap2_seq else { continue; };

            // assign ambiguous reads to a haplotype
            let ambiguous_read_is: Vec<ReadIndex> = read_haplotypes.iter().enumerate()
                .filter_map(|(read_i, hap)|{
                    if *hap == 0 { Some(read_i) } else { None }
                }).collect();
            if ambiguous_read_is.len() > 0 {
                ambiguous_read_is.iter().for_each(|read_i| {
                    read_haplotypes[*read_i] = self.assign_ambiguous_to_haplotype(
                        &reads[*read_i], &hap1_seq, &hap2_seq
                    );
                });    
                hap1_read_is = read_haplotypes.iter().enumerate()
                    .filter_map(|(read_i, hap)|{
                        if *hap == 1 { Some(read_i) } else { None }
                    }).collect();
                hap2_read_is = read_haplotypes.iter().enumerate()
                    .filter_map(|(read_i, hap)|{
                        if *hap == 2 { Some(read_i) } else { None }
                    }).collect(); 
                // TODO: regenerate consensus including newly assigned ambiguous reads?
            }

            // align each read to its haplotype consensus to call subclonal variants
            self.align_to_haplotype_consensus(
                &mut op_val, &mut ref_pos0_map, &mut hap_vs_ref, 
                haplotype_consensuses, &mut frag_vars, &mut encoding,
                re_fragment, reads, 
                hap1_read_is, hap1_seq, Haplotype::Haplotype1
            );
            self.align_to_haplotype_consensus(
                &mut op_val, &mut ref_pos0_map, &mut hap_vs_ref, 
                haplotype_consensuses, &mut frag_vars, &mut encoding, 
                re_fragment, reads, 
                hap2_read_is, hap2_seq, Haplotype::Haplotype2
            );
        }
    }

    /// Use unambiguous haplotype reads to build a haplotype consensus sequence.
    fn build_haplotype_consensus(
        &mut self,
        seq0_bases: &mut Vec<String>,
        cs_map:     &mut Vec<FxHashMap<String, u8>>,
        op_val:     &mut String,
        reads:      &[ReadInstance],
        read_is:    &[ReadIndex],
    ) -> Option<UppercaseACGTN> {

        // initialize the first forward strand sequence as the abitrary target
        // remembering that nearly all reads are forward strand after basecalling
        let read0_i = read_is.iter().position(|i| !reads[*i].is_reverse);
        let Some(read0_i) = read0_i else { return None; }; // never expected to fail
        let seq0_bytes = &reads[read0_i].seq_bytes;
        let mm2 = self.minimap2.clone()
            .with_seq(seq0_bytes)
            .expect("Failed to initialize minimap2 in build_haplotype_consensus()");
        seq0_bases.clear();
        seq0_bases.extend(seq0_bytes.iter()
            .map(|&b| match b {
                b'A' => "A".to_string(),
                b'C' => "C".to_string(),
                b'G' => "G".to_string(),
                b'T' => "T".to_string(),
                _    => "N".to_string(),
            })
        );
        cs_map.clear();
        cs_map.extend(seq0_bases.iter()
            .map(|b| {
                let mut m = FxHashMap::default();
                m.insert(b.clone(), 1);
                m
            })
        );

        // align each remaining read to target to count alternative bases
        read_is.iter().for_each(|read_i| {
            if *read_i == read0_i { return; } // skip the read in use as target
            let mappings = mm2.map(
                &reads[*read_i].seq_bytes, 
                true, 
                false, 
                None, 
                Some(&MM_F_NO_PRINT_2ND), 
                None
            ).expect("Minimap2 failed in build_haplotype_consensus()");
            for mapping in mappings {
                let Some(aln) = mapping.alignment else { continue; };
                let Some(cs) = aln.cs else { continue; };
                let mut tgt_pos0 = mapping.target_start as usize;
                let mut chars = cs.chars();
                let mut op = chars.next().unwrap();
                op_val.clear();
                while let Some(char) = chars.next() {
                    if char.is_alphanumeric() {
                        op_val.push(char);
                    } else {
                        match op {
                            ':' => { // :[0-9]+   Identical sequence length
                                (0..op_val.parse::<usize>().unwrap()).for_each(|_| {
                                    *cs_map[tgt_pos0].get_mut(&seq0_bases[tgt_pos0]).unwrap() += 1;
                                    tgt_pos0 += 1;
                                });
                            },
                            '*' => { // *[acgtn][acgtn]   Substitution: target to query
                                cs_map[tgt_pos0]
                                    .entry(op_val[1..=1].to_ascii_uppercase())
                                    .and_modify(|n| *n += 1)
                                    .or_insert(1);
                                tgt_pos0 += 1;
                            },
                            '+' => { // +[acgtn]+   Insertion to the target
                                let alt = op_val.to_ascii_uppercase();
                                cs_map[tgt_pos0 - 1]
                                    .entry(format!("{}{}", seq0_bases[tgt_pos0 - 1], alt))
                                    .and_modify(|n| *n += 1)
                                    .or_insert(1);
                            },
                            '-' => { // -[acgtn]+   Deletion from the target
                                (0..op_val.len()).for_each(|_| {
                                    cs_map[tgt_pos0]
                                        .entry("-".to_string())
                                        .and_modify(|n| *n += 1)
                                        .or_insert(1);
                                    tgt_pos0 += 1;
                                });
                            },
                            _ => panic!("Unexpected CS tag operation: {}", op),
                        }
                        op = char;
                        op_val.clear();
                    }
                }
                if let Ok(len) = op_val.parse::<usize>() {
                    (0..len).for_each(|_| {
                        *cs_map[tgt_pos0].get_mut(&seq0_bases[tgt_pos0]).unwrap() += 1;
                        tgt_pos0 += 1;
                    });
                }
            }
        });

        // scan the matrix by base to establish the haplotype consensus
        let n_seq0_bases = cs_map.len();
        let mut consensus = String::with_capacity(n_seq0_bases + 100);
        for tgt_pos0 in 0..n_seq0_bases {
            let bases = cs_map[tgt_pos0].iter()
                .max_by(|a, b| a.1.cmp(&b.1))
                .map(|(bases, _)| bases)
                .unwrap();
            if bases != "-" {
                consensus.push_str(bases); // usually one but could be multiple bases at an insertion
            }
        }
        Some(consensus)
    }

    /// Align one ambiguous read to each haplotype consensus to resolve its 
    /// haplotype.
    fn assign_ambiguous_to_haplotype(
        &mut self,
        read:     &ReadInstance,
        hap1_seq: &UppercaseACGTN,
        hap2_seq: &UppercaseACGTN,
    ) -> u8 {
        let mm2 = self.minimap2.clone()
            .with_seq(&read.seq_bytes)
            .expect("Failed to initialize minimap2 in assign_ambiguous_to_haplotype()");
        let m1 = mm2.map(
            hap1_seq.as_bytes(), 
            true, // cs not used, but need full score calculation
            false, 
            None, 
            Some(&MM_F_NO_PRINT_2ND), 
            None
        ).expect("Minimap2 failed in assign_ambiguous_to_haplotype()");
        let m2 = mm2.map(
            hap2_seq.as_bytes(), 
            true, 
            false, 
            None, 
            Some(&MM_F_NO_PRINT_2ND), 
            None
        ).expect("Minimap2 failed in assign_ambiguous_to_haplotype()");
        let score1 = m1[0].alignment.as_ref()
            .map_or(0, |a| {
                a.alignment_score.unwrap_or(0)
            });
        let score2 = m2[0].alignment.as_ref()
            .map_or(0, |a| {
                a.alignment_score.unwrap_or(0)
            });
        if score1 > score2 { 1 } else { 2 }
    }

    /// Align all haplotype read sequences to their consensus to call subclonal
    /// variant relative to that consensus.
    fn align_to_haplotype_consensus(
        &mut self,
        op_val:       &mut String,
        ref_pos0_map: &mut Vec<ChromPos0>,
        hap_vs_ref:   &mut String,
        haplotype_consensuses: &mut HaplotypeConsensuses,
        frag_vars:   &mut FragmentVariants,
        encoding:    &mut AlignmentEncoding,
        re_fragment: &ReFragment,
        reads:       &[ReadInstance],
        read_is:     Vec<ReadIndex>,
        hap_seq:     UppercaseACGTN,
        haplotype:   Haplotype,
    ) {
        let n_hap_bases = hap_seq.len();
        let mm2 = self.minimap2.clone()
            .with_seq(hap_seq.as_bytes())
            .expect("Failed to initialize minimap2 in align_to_haplotype_consensus()");

        // align the fragment reference span to the haplotype consensus, i.e., ref_on_hap
        let (ref_seq, _) = haplotype_consensuses.get(re_fragment, Haplotype::Unspecified);
        let ref_on_hap = &mm2.map(
            ref_seq.as_bytes(), 
            true, 
            false, 
            None, 
            Some(&MM_F_NO_PRINT_2ND), 
            None
        ).expect("Minimap2 failed to align ref_on_hap in align_to_haplotype_consensus()")[0];
        let Some(aln) = &ref_on_hap.alignment else { return; };
        let Some(cs) = &aln.cs else { return; };

        // create a map of reference positions per each haplotype consensus position
        // and a hap_vs_ref encoding to retain a memory of clonal variants relative to reference
        ref_pos0_map.clear();
        hap_vs_ref.clear();
        if ref_on_hap.target_start > 0 {
            let count = ref_on_hap.target_start as usize;
            ref_pos0_map.extend(repeat_n(re_fragment.start0, count));
            hap_vs_ref.extend(repeat_n("+", count));
        }
        let mut ref_pos0 = re_fragment.start0 + ref_on_hap.query_start as u32;
        let mut chars = cs.chars();
        let mut op = chars.next().unwrap();
        op_val.clear();
        while let Some(char) = chars.next() {
            if char.is_alphanumeric() {
                op_val.push(char);
            } else {
                match op {
                    ':' => { // :[0-9]+   Identical sequence length
                        let len = op_val.parse::<usize>().unwrap();
                        (0..len).for_each(|_| {
                            ref_pos0_map.push(ref_pos0);
                            ref_pos0 += 1;
                        });
                        hap_vs_ref.push_str(&format!("={}", len));
                    },
                    '*' => { // *[acgtn][acgtn]   Substitution: target to query
                        ref_pos0_map.push(ref_pos0);
                        ref_pos0 += 1;
                        hap_vs_ref.push_str(&op_val[0..=0].to_ascii_uppercase());
                    },
                    '+' => { // +[acgtn]+   Insertion to the target
                        hap_vs_ref.pop();
                        hap_vs_ref.push('-');
                        ref_pos0 += op_val.len() as u32;
                    },
                    '-' => { // -[acgtn]+   Deletion from the target
                        (0..op_val.len()).for_each(|_| {
                            ref_pos0_map.push(ref_pos0 - 1);
                            hap_vs_ref.push('+');
                        });
                    },
                    _ => panic!("Unexpected CS tag operation: {}", op),
                }
                op = char;
                op_val.clear();
            }
        }
        if let Ok(len) = op_val.parse::<usize>(){
            (0..len).for_each(|_| {
                ref_pos0_map.push(ref_pos0);
                ref_pos0 += 1;
            }); 
            hap_vs_ref.push_str(&format!("={}", len));
        }
        if ref_on_hap.target_end < n_hap_bases as i32 - 1 {
            let count = n_hap_bases - 1 - ref_on_hap.target_end as usize;
            ref_pos0_map.extend(repeat_n(re_fragment.end1 - 1, count));
            hap_vs_ref.extend(repeat_n("+", count));
        }

        // align each haplotype read to its consensus
        let n_reads = read_is.len();
        frag_vars.reset(n_reads);
        for read_j in 0..n_reads {
            let read = &reads[read_is[read_j]];
            let read_on_hap = &mm2.map(
                &read.seq_bytes, 
                true, 
                false, 
                None, 
                Some(&MM_F_NO_PRINT_2ND), 
                None
            ).expect("Minimap2 failed to align read_on_hap in align_to_haplotype_consensus()")[0];

            // extract each read's subclonal variants and commit its haplotype encoding
            if let Some(aln) = &read_on_hap.alignment  {
                if let Some(cs) = &aln.cs  {
                    encoding.prepare_read_on_hap(re_fragment, read, read_on_hap.target_start as usize);
                    self.process_cs_tag(
                        haplotype, 
                        Some(read_on_hap), Some(cs), Some(&ref_pos0_map),
                        frag_vars, encoding, 
                        read_j, read
                    );
                } else {
                    // unexpected alignment failure, mask the entire encoding
                    encoding.prepare_read_on_hap(re_fragment, read, n_hap_bases);
                }
            } else {
                encoding.prepare_read_on_hap(re_fragment, read, n_hap_bases);
            }
            self.reads_on_haplotype.insert(
                re_fragment, 
                haplotype,
                encoding.clone()
            ); 
        }

        // call subclonal variants aggregated over all haplotype reads (might be none)
        for hap_var in frag_vars.variant_map.keys(){
            let vmap =  &frag_vars.variant_map[&hap_var];
            let read_js: Vec<ReadIndex> = vmap.read_map.iter()
                .enumerate()
                .filter_map(|(read_j, r)|{
                    if r.has_var() { 

                        Some(read_j) 
                    } else { None }
                }).collect();
            let max_avg_qual = read_js.iter()
                .map(|read_j| {
                    let avg_qual= vmap.read_map[*read_j].avg_qual;
                    self.variant_reads_tally.add_subclonal_variant(
                        &reads[read_is[*read_j]], re_fragment, &haplotype, 
                        hap_var, avg_qual
                    );
                    avg_qual
                })
                .max()
                .unwrap_or_default();
            self.variant_tally.add_subclonal(
                &hap_var, reads, &read_is, &read_js, 
                max_avg_qual
            );
        }

        // cache the haplotype consensus for printing
        haplotype_consensuses.insert(
            re_fragment, haplotype, 
            hap_seq, Some(hap_vs_ref.clone())
        );
    }
}

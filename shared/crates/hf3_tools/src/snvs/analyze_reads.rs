//! Support for building haplotype consensuses from RE fragments and using them
//! to call rare variants without reference alignment bias.
//! 
//! Clonal variants are called from recurring variants seen in prior minimap2
//! alignments against reference. Subclonal variants are called from new
//! Smith-Waterman alignments against haplotype consensus(es) generated here.

// dependencies
use std::iter::repeat_n;
use rustc_hash::FxHashMap;
use genomex::sequence::AlignmentStatus;
use super::*;

// constants
const MIN_HAPLOTYPE_READS:    usize = 3;  // a haplotype must have >=3 matching reads
const MIN_HETEROZYGOUS_READS: usize = 8;  // allow fewer total reads when heterozygous
const MIN_HOMOZYGOUS_READS:   usize = 15; // require more homozygous reads to minimize false homozygosity

impl SnvChromWorker {

    /// Build homozygous or heterozygous consensuses for RE fragment haplotypes,
    /// and use them to call clonal and subclonal variants.
    pub fn analyze_reads(
        &mut self,
        tool: &SnvAnalysisTool,
        fragment_reads: FragmentReads,
        haplotype_consensuses: &mut HaplotypeConsensuses,
    ){
        // allocate required objects
        let mut frag_vars = FragmentVariants::new();
        let mut encoding = AlignmentEncoding::new();

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
                    let read_is: Vec<ReadIndex> = vmap.read_map.iter().enumerate()
                        .filter_map(|(read_i, b)|{
                            if *b { Some(read_i) } else { None }
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

            // require more reads to continue with a homozygous RE fragment
            let n_heterozygous_variants = frag_vars.variant_map.len();
            if n_heterozygous_variants == 0 {
                if n_reads < MIN_HOMOZYGOUS_READS { continue; }
                let hap3_read_is: Vec<ReadIndex> = (0..n_reads).map(|i| i). collect();
                let hap3_seq = self.build_haplotype_consensus(
                    &reads, &hap3_read_is
                );
                self.align_to_haplotype_consensus(
                    haplotype_consensuses, &mut frag_vars, &mut encoding, 
                    re_fragment, reads, 
                    hap3_read_is, hap3_seq, Haplotype::Homozygous
                );
                continue;
            }

            // normalize non-reference variants to consistent haplotype numbers
            // given that variant SNPs may be on opposing haplotypes
            // reads that cannot be unambiguously assigned are haplotype 0 at this point
            let het_vars: Vec<_> = frag_vars.variant_map.keys().cloned().collect();
            let index_matches = frag_vars.variant_map[&het_vars[0]].read_map.clone();
            let mut read_haplotypes: Vec<u8> = if n_heterozygous_variants > 1 { 
                // re-orient variant read maps to index_matches, i.e., to the first variant
                for het_var in het_vars[1..n_heterozygous_variants].iter(){
                    let vmap =  frag_vars.variant_map.get_mut(&het_var).unwrap();
                    let match_score = vmap.read_map.iter().zip(&index_matches).map(|b|{
                        if *b.0 == *b.1 { 1 } else { -1 }
                    }).sum::<i32>();
                    if match_score < 0 {
                        for b in vmap.read_map.iter_mut() { *b = !*b; }
                    }
                }
                // assign each unambiguous read to haplotype 1 or 2; 0 if ambiguous
                (0..n_reads).map(|read_i|{
                    let rmap: Vec<bool> = het_vars.iter().map(|het_var|{
                        let vmap =  frag_vars.variant_map.get(&het_var).unwrap();
                        vmap.read_map[read_i]
                    }).collect();
                    if rmap.windows(2).all(|w| w[0] == w[1]) {
                        rmap[0] as u8 + 1
                    } else {
                        0 // a read wasn't always assigned to the same haplotype
                    }
                }).collect()
            } else {
                // encoding assigns the variant allele as haplotype2 for single SNPs
                index_matches.into_iter().map(|b| b as u8 + 1).collect()
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
            let hap1_seq = self.build_haplotype_consensus(&reads, &hap1_read_is);
            let hap2_seq = self.build_haplotype_consensus(&reads, &hap2_read_is);

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
                haplotype_consensuses, &mut frag_vars, &mut encoding,
                re_fragment, reads, 
                hap1_read_is, hap1_seq, Haplotype::Haplotype1
            );
            self.align_to_haplotype_consensus(
                haplotype_consensuses, &mut frag_vars, &mut encoding, 
                re_fragment, reads, 
                hap2_read_is, hap2_seq, Haplotype::Haplotype2
            );
        }
    }

    /// Use unambiguous haplotype reads to build a haplotype consensus sequence.
    fn build_haplotype_consensus(
        &mut self,
        reads:   &[ReadInstance],
        read_is: &[ReadIndex],
    ) -> UppercaseACGTN {

        // initialize the first sequence as the abitrary target
        let seq0 = reads[read_is[0]].get_top_strand_seq();
        let n_bases = seq0.len();
        let n_reads = read_is.len();
        let mut map: Vec<Vec<UppercaseACGTN>> = Vec::with_capacity(n_reads);
        map.push(seq0.split("").map(|c| c.to_string()).collect());

        // align each remaining read to target to fill the matrix
        for read_j in 1..n_reads{
            let aln = self.aligner.align(
                &reads[read_is[read_j]].get_top_strand_seq(), 
                &seq0, 
                None, 
                false,
            );
            if aln.status == AlignmentStatus::AlignmentFound {
                if aln.tgt_start0 > 0 || aln.tgt_end0 < n_bases - 1 {
                    let mut qry_on_tgt = Vec::with_capacity(n_bases);
                    if aln.tgt_start0 > 0 {
                        qry_on_tgt.extend(repeat_n("N".to_string(), aln.tgt_start0));
                    }
                    qry_on_tgt.extend(aln.qry_on_tgt.into_iter());
                    if aln.tgt_end0 < n_bases - 1 {
                        qry_on_tgt.extend(repeat_n("N".to_string(), n_bases - 1 - aln.tgt_end0));
                    }
                    map.push(qry_on_tgt);   
                } else {
                    map.push(aln.qry_on_tgt);   
                }
            }
        }

        // scan the matrix by base to establish the haplotype consensus
        let n_reads = map.len(); // just in case one failed to align, not expected
        (0..n_bases)
            .map(|base_i|{

                // determine the most frequently observed value at each index base position
                // could be a single base, multiple bases, or "-" for deletion
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

            // remove consensus deletion ops when index target had an insertion
            .filter_map(|base|{ 
                if base == "-" { None } else { Some(base.to_owned()) }
            })
            .collect::<Vec<String>>()

            // concatenate to the final consensus sequence
            .join("")
    }

    /// Align one ambiguous read to each haplotype consensus to resolve its 
    /// haplotype.
    fn assign_ambiguous_to_haplotype(
        &mut self,
        read: &ReadInstance,
        hap1_seq: &UppercaseACGTN,
        hap2_seq: &UppercaseACGTN,
    ) -> u8 {
        let seq = read.get_top_strand_seq();
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

    /// Align all haplotype read sequences to their consensus to call subclonal
    /// variant relative to that consensus.
    fn align_to_haplotype_consensus(
        &mut self,
        haplotype_consensuses: &mut HaplotypeConsensuses,
        frag_vars:   &mut FragmentVariants,
        encoding:    &mut AlignmentEncoding,
        re_fragment: &ReFragment,
        reads:       &[ReadInstance],
        read_is:     Vec<ReadIndex>,
        hap_seq:     UppercaseACGTN,
        haplotype:   Haplotype,
    ) {
        // align the fragment reference span to the haplotype consensus 
        // to enable mapping the simple repeat map onto the haplotype consensus
        let aln = self.aligner.align(
            haplotype_consensuses.get(re_fragment, Haplotype::Unspecified), 
            &hap_seq, 
            None, 
            false,
        );
        if aln.status != AlignmentStatus::AlignmentFound { return; }
        let n_bases = hap_seq.len();
        let mut ref_pos0_map: Vec<ChromPos0> = Vec::with_capacity(n_bases);
        if aln.tgt_start0 > 0 {
            ref_pos0_map.extend(repeat_n(re_fragment.start0, aln.tgt_start0));
        }
        let mut ref_pos0 = re_fragment.start0 + aln.qry_start0 as u32;
        for val in aln.qry_on_tgt {
            if val == "-" {
                ref_pos0_map.push(ref_pos0 - 1);
            } else if val.len() > 1 {
                ref_pos0 += val.len() as u32 - 1;
                ref_pos0_map.push(ref_pos0);
                ref_pos0 += 1;
            } else {
                ref_pos0_map.push(ref_pos0);
                ref_pos0 += 1;
            }
        }
        if aln.tgt_end0 < n_bases - 1 {
            ref_pos0_map.extend(repeat_n(re_fragment.end1 - 1, n_bases - 1 - aln.tgt_end0));
        }

        // align each haplotype read to its consensus
        let n_reads = read_is.len();
        frag_vars.reset(n_reads);
        for read_j in 0..n_reads {
            let read = &reads[read_is[read_j]];
            let aln = self.aligner.align(
                &read.get_top_strand_seq(), 
                &hap_seq, 
                None, 
                false,
            );
            if aln.status == AlignmentStatus::AlignmentFound {

                // extract read subclonal variants and commit its haplotype encoding
                encoding.prepare_read_on_hap(re_fragment, read, aln.tgt_start0);
                let cs = SnvChromWorker::get_cs_tag(&aln, &hap_seq);
                self.process_cs_tag(
                    haplotype, 
                    Some(aln), Some(cs), Some(&ref_pos0_map),
                    frag_vars, encoding, 
                    read_j, read
                );
                self.reads_on_haplotype.insert(
                    re_fragment, 
                    haplotype,
                    encoding.clone()
                );
            }
        }

        // call subclonal variants aggregated over all haplotype reads
        for hap_var in frag_vars.variant_map.keys(){
            let vmap =  &frag_vars.variant_map[&hap_var];
            let qmap = &frag_vars.qual_map;
            let read_js: Vec<ReadIndex> = vmap.read_map.iter().enumerate()
                .filter_map(|(read_j, b)|{
                    if *b { Some(read_j) } else { None }
                }).collect();
            let max_avg_qual = read_js.iter()
                .map(|read_j| qmap[&(hap_var.clone(), *read_j)])
                .max()
                .unwrap_or_default();
            self.variant_tally.add_subclonal(
                &hap_var, reads, &read_is, &read_js, 
                max_avg_qual
            );
        }

        // cache the haplotype consensus for printing
        haplotype_consensuses.insert(re_fragment, haplotype, hap_seq);
    }
}

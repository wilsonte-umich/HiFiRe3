//  action:
//    junctions are rejected if:
//        the flanking alignments are consistent with a foldback inversion 
//        they are from ONT and consistent with a follow-on event with low-quality inserted bases
//        an adapter is present in junction inserted bases of sufficient length
//    return a junction failure flag
//  note: additional chimeric splitting occurs after final RE site assessment downstream

// dependencies
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use genomex::sam::{SamRecord, flag, junction::Junction};
use genomex::sequence::{Aligner, rc_acgt_str};
use crate::junctions::JxnFailureFlag;

// constants
pub_key_constants!{
    // from environment variables
    INSERTION_WINDOW_SIZE
    MIN_INSERTION_WINDOW_QUAL
    INSERTION_ADAPTER_SEQUENCE
    MIN_ADAPTER_LENGTH
    MAX_ADAPTER_LENGTH
    MIN_ADAPTER_SCORE
    // counters
    ADAPTER_SCORES
    N_CHIMERIC
}
const PHRED_OFFSET: usize = 33;

pub struct ChimeraSplitter {
    is_ont: bool,
    insertion_window_size: isize,
    min_sum_ins_qual: usize,
    has_adapters: bool,
    adapter_core: String,
    adapter_core_rc: String,
    min_adapter_length: isize,
    max_adapter_length: isize,
    min_adapter_score: i32,
    max_possible_score: i32,
    aligner: Aligner,
}
impl ChimeraSplitter {
    /* ---------------------------------------------------------------------------
    initialize
    ---------------------------------------------------------------------------- */
    /// Initialize a junction chimera splitter.
    pub fn new(w: &mut Workflow) -> ChimeraSplitter {
        w.cfg.set_usize_env(&[INSERTION_WINDOW_SIZE, MIN_INSERTION_WINDOW_QUAL, 
                                    MIN_ADAPTER_LENGTH, MAX_ADAPTER_LENGTH, MIN_ADAPTER_SCORE]);
        w.cfg.set_string_env(&[INSERTION_ADAPTER_SEQUENCE]);
        let insertion_window_size = *w.cfg.get_usize(INSERTION_WINDOW_SIZE);
        let min_insertion_window_qual = *w.cfg.get_usize(MIN_INSERTION_WINDOW_QUAL);
        let adapter_core = w.cfg.get_string(INSERTION_ADAPTER_SEQUENCE); 
        let has_adapters = !adapter_core.is_empty();
        let min_adapter_length = *w.cfg.get_usize(MIN_ADAPTER_LENGTH);
        let max_adapter_length = *w.cfg.get_usize(MAX_ADAPTER_LENGTH);
        let aligner = Aligner::new(
            adapter_core.len().max(min_adapter_length), 
            max_adapter_length
        ).suppress_alignment_map();
        w.ctrs.add_counters(&[
            (N_CHIMERIC, "chimeric reads identified by junction analysis"),
        ]);
        if has_adapters {
            w.ctrs.add_indexed_counters(&[
                (ADAPTER_SCORES, "sw_score", adapter_core.len(), "Smith-Waterman alignment scores for adapter vs. inserted bases"),
            ]);
        }
        ChimeraSplitter{
            is_ont: *w.cfg.get_bool(super::IS_ONT), // must be set upstream
            insertion_window_size: insertion_window_size as isize,
            min_sum_ins_qual: (min_insertion_window_qual + PHRED_OFFSET) * insertion_window_size,
            has_adapters,
            adapter_core:      adapter_core.to_string(),       // fused to 5' genomic ends; for ONT ligation kit, last T matches the one-base A-tail
            adapter_core_rc:   rc_acgt_str(&adapter_core), // fused to 3' genomic ends 
            min_adapter_length: min_adapter_length as isize,
            max_adapter_length: max_adapter_length as isize,
            min_adapter_score: *w.cfg.get_usize(MIN_ADAPTER_SCORE) as i32,
            max_possible_score: adapter_core.len() as i32,
            aligner: aligner,
        }
    }

    /// check junction quality from pairs of alignments
    pub fn get_jxn_failure_flag(
        &mut self, 
        w: &mut Workflow, 
        alns: &[SamRecord],
        aln5_i: usize, 
        aln3_i: usize,
        jxn: &Junction,
    ) -> JxnFailureFlag {

        // reject junctions consistent with foldback inversions
        // do not increment chimeric count for foldbacks, they are intra, not inter-molecular
        if Self::is_foldback(&alns[aln5_i], &alns[aln3_i]) {
            w.ctrs.increment_keyed(super::N_JXNS_BY_REASON, super::JXN_FAIL_FOLDBACK_INV);
            return JxnFailureFlag::FoldbackInv;
        }

        // reject ONT junctions consistent with follow-on events
        // let (clip1, jxn_ins_size, is_follow_on) = self.is_ont_follow_on(aln5, aln3);
        if self.is_ont && self.is_ont_follow_on(jxn) {
            w.ctrs.increment(N_CHIMERIC);
            w.ctrs.increment_keyed(super::N_JXNS_BY_REASON, super::JXN_FAIL_ONT_FOLLOW_ON);
            return JxnFailureFlag::OntFollowOn;
        }

        // reject junctions that contain adapters
        if self.jxn_has_adapters(w, jxn){
            w.ctrs.increment(N_CHIMERIC);
            w.ctrs.increment_keyed(super::N_JXNS_BY_REASON, super::JXN_FAIL_HAS_ADAPTER);
            return JxnFailureFlag::HasAdapter;
        }

        // junction passed all criteria
        w.ctrs.increment_keyed(super::N_JXNS_BY_REASON, super::JXN_FAIL_NONE);
        JxnFailureFlag::None

    }

    // /* ---------------------------------------------------------------------------
    // chimera filter functions
    // ---------------------------------------------------------------------------- */
    /// Check whether a sequence is consistent with an ONT duplex read as a single 
    /// foldback inversion. It is most sensitive and acceptable to reject any 
    /// inversion junction with reverse-complement overlap between its flanking alignments.
    /// 
    /// UPDATE: This check now performed for all library types, not just ONT.
    //// ```
    /// ----->
    ///      | inversion junction
    /// <-----
    /// ```
    fn is_foldback(aln5: &SamRecord, aln3: &SamRecord) -> bool {
        if aln5.rname != aln3.rname { return false; } // translocation
        if aln5.check_flag_any(flag::REVERSE) == 
           aln3.check_flag_any(flag::REVERSE) { return false; } // deletion or duplication
        if aln5.pos1 > aln3.get_end1() || 
           aln3.pos1 > aln5.get_end1(){ return false; } // no overlap
        true
    }

    /// Determine if an ONT junction has a very low quality insertion span,
    /// identifying it as a two-insert follow-on event.
    fn is_ont_follow_on(&self, jxn: &Junction) -> bool {

        // stop if insertion of insufficient size to warrant adapter/follow-on detection
        if jxn.offset <  self.min_adapter_length ||
           jxn.offset <= self.insertion_window_size {
            return false;
        }

        // examine windows of bases to find no-base signal stretches
        // observed behavior is read1 ... (3' adapter) ... no-base signals/bases ... 5' adapter ... read2
        let quals = jxn.jxn_qual.as_bytes();
        let mut sum_ins_qual: usize = quals[..self.insertion_window_size as usize].iter().map(|&b| b as usize).sum();
        if sum_ins_qual < self.min_sum_ins_qual {
            return true; // the first window is low quality
        }
        for i in 1..=(jxn.offset - self.insertion_window_size) {
            sum_ins_qual -= quals[i as usize - 1] as usize;
            sum_ins_qual += quals[(i + self.insertion_window_size - 1) as usize] as usize;
            if sum_ins_qual < self.min_sum_ins_qual {
                return true; // a later window is low quality
            }
        }
        false
    }

    /// Determine if a junction insertion has adapters, identifying it as a two-insert event.
    fn jxn_has_adapters(&mut self, w: &mut Workflow,jxn: &Junction,) -> bool {

        // stop if no adapter sequences was provided for the library type
        // or if insertion of insufficient size to warrant adapter detection
        if !self.has_adapters ||
           jxn.offset < self.min_adapter_length ||
           jxn.offset > self.max_adapter_length {
            return false;
        }

        // next examine the inserted bases for adapters using Smith-Waterman on both strands
        let sw = self.aligner.align(&self.adapter_core, &jxn.jxn_seq, None, true);
        w.ctrs.increment_indexed(ADAPTER_SCORES, sw.score.max(0).min(self.max_possible_score) as usize);
        if sw.score >= self.min_adapter_score {
            return true;
        }
        let sw = self.aligner.align(&self.adapter_core_rc, &jxn.jxn_seq, None, true);
        w.ctrs.increment_indexed(ADAPTER_SCORES, sw.score.max(0).min(self.max_possible_score) as usize);
        if sw.score >= self.min_adapter_score {
            return true;
        }
        false
    }

    // /* ---------------------------------------------------------------------------
    // export likely SV artifact reads for exploring artifact mechanisms
    // ---------------------------------------------------------------------------- */
    // /// export likely SV artifact reads to a separate file
    // /// only on-target reads with exactly two alignments reach this sub
    // sub recordVariantSizes {
    //     my ($sizeBin, @alns) = @_;

    //     // do not record foldback inversions as intermolecular chimeras
    //     $isFoldback and return;
    //     my $stemLengthBin = int(min(
    //         $alns[0][STEM5_LENGTH],
    //         $alns[1][STEM3_LENGTH] ? $alns[1][STEM3_LENGTH] : 1e9
    //     ) / $sizePlotBinSize) * $sizePlotBinSize;

    //     // record insert sizes for all single-junction SV reads stratified by chimericity
    //     // non-chimeric here may contain deletion and other true SVs
    //     if($isChimeric){
    //         $insertSizeCounts{$READ_LEN_CHIMERIC}{$sizeBin}++;
    //         $stemLengthCounts{$READ_LEN_CHIMERIC}{$stemLengthBin}++;
    //     } else {
    //         $insertSizeCounts{$READ_LEN_NON_CHIMERIC}{$sizeBin}++;
    //         $stemLengthCounts{$READ_LEN_NON_CHIMERIC}{$stemLengthBin}++;
    //     }

    //     // only print single-junction interchromosomal chimeras as presumptive artifacts
    //     // we do not yet know if these are single-molecule or would be confirmed downstream
    //     // but most interchromosomal, i.e., translocation molecules are likely artifacts
    //     $alns[0][S_RNAME] eq $alns[1][S_RNAME] and return;

    //     // stratify interchromosomal molecules as inter- or intra-genomic, when applicable
    //     my $genomicity = $IS_COMPOSITE_GENOME ?
    //         isInterGenomic($alns[0][S_RNAME], $alns[1][S_RNAME]) :
    //         INTRAGENOMIC;

    //     // reads were stratified above as chimeric or not
    //     // "chimeric" means a single junction that failed breakpointMatchesSite, isOntFollowOn, or junctionHasAdapters
    //     // "not chimeric" means it passed those chimeric tests, but it may still be chimeric by some other mechanism
    //     my $chimericity = $isChimeric ? CHIMERIC : NOT_CHIMERIC;
    //     my $type = $genomicity.CHIMERIC_DELIMITER.$chimericity;

    //     // increment size counts for translocation presumptive artifacts
    //     $insertSizeCounts{$type}{$sizeBin}++;
    //     $stemLengthCounts{$type}{$stemLengthBin}++;
    // }
    // sub isInterGenomic {
    //     my ($chrom1, $chrom2) = @_;
    //     my ($genome1) = ($chrom1 =~ m/.+_(.+)/); // e.g., chr1_(hs1)
    //     my ($genome2) = ($chrom2 =~ m/.+_(.+)/);
    //     $genome1 eq $genome2 ? INTRAGENOMIC : INTERGENOMIC;
    // }

    // /// print summary information
    // sub printChimeraSummary {

    //     // print and save insert size distributions
    //     open my $sizesH, '>', "$FILTERED_INSERT_SIZES_FILE" or throwError("could not open insert sizes file: $!");
    //     open my $stemsH, '>', "$FILTERED_STEM_LENGTHS_FILE" or throwError("could not open stem lengths file: $!");
    //     foreach my $type(
    //         READ_LEN_PROPER, PROJ_LEN_PROPER, READ_LEN_CHIMERIC, READ_LEN_NON_CHIMERIC, 
    //         @translocationTypes
    //     ){
    //         if($insertSizeCounts{$type}){
    //             my @insertSizes = sort { $a <=> $b } keys %{$insertSizeCounts{$type}};
    //             print $sizesH join("\n", map { 
    //                 join("\t", 
    //                     $type,
    //                     $_, 
    //                     $insertSizeCounts{$type}{$_}
    //                 )
    //             } @insertSizes), "\n";
    //         }
    //         if($stemLengthCounts{$type}){
    //             my @stemLengths = sort { $a <=> $b } keys %{$stemLengthCounts{$type}};
    //             print $stemsH join("\n", map { 
    //                 join("\t", 
    //                     $type,
    //                     $_, 
    //                     $stemLengthCounts{$type}{$_}
    //                 )
    //             } @stemLengths), "\n";
    //         }
    //     }
    //     close $sizesH;
    //     close $stemsH;

    //     // print site distance distribution to log
    //     if($REJECTING_JUNCTION_RE_SITES){
    //         print STDERR "\njunction node-to-site distances\n";
    //         foreach my $dist(0..MAX_SITE_DISTANCE){
    //             print STDERR join("\t", $dist, $siteDistances[$dist] || 0), "\n";
    //         }
    //     }

    //     // print adapter SW score distribution to log
    //     print STDERR "\nadapter SW scores\n";
    //     foreach my $score(0..100){
    //         print STDERR join("\t", $score, $adapterScores[$score] || 0), "\n";
    //     }
    // }

}

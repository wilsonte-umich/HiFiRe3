//! action:
//!   report insert sizes of various types

// dependencies
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use genomex::sam::{SamRecord, flag};
use crate::junctions::{JxnFailureFlag, is_chimeric_jxn};
use super::Tool;

// constants
pub_key_constants!{
    // from environment variables
    IS_COMPOSITE_GENOME
    FILTERED_INSERT_SIZES_FILE
    FILTERED_STEM_LENGTHS_FILE
    // counters
    INSERT_SIZES_CHIMERIC
    INSERT_SIZES_NON_CHIMERIC
    STEM_LENGTHS_CHIMERIC
    STEM_LENGTHS_NON_CHIMERIC
}
const READ_LEN_PROPER:       &str = "nonSV.actual";
const PROJ_LEN_PROPER:       &str = "nonSV.projected";
const READ_LEN_CHIMERIC:     &str = "SV.chimeric";   // cannot assess projected size of structural variant reads
const READ_LEN_NON_CHIMERIC: &str = "SV.passed"; // these include all SV reads, not just translocations
// -------------
const INTRAGENOMIC: &str = "intra"; // for assembling types of translocation, i.e., presumed artifacts
const INTERGENOMIC: &str = "inter";
const CHIMERIC:     &str = "chimeric";
const NOT_CHIMERIC: &str = "passed";
const CHIMERIC_DELIMITER: &str = ".";

    // N_BY_STEM_LENGTH
    // MIN_SELECTED_SIZE
    // SELECTED_SIZE_CV
    // MIN_ALLOWED_SIZE
    // HAS_BASE_ACCURACY
    // PLATFORM_MAX_INSERT_SIZE


/// InsertSizer helps analyze insert size and junction stem length distributions.
pub struct InsertSizer {
    is_composite_genome: bool,
    translocation_types: Vec<String>,
    insert_size: usize,
    stem_length: usize,
}
impl InsertSizer {
    /* ---------------------------------------------------------------------------
    initialize
    ---------------------------------------------------------------------------- */
    /// Initialize an insert size recorder.
    pub fn new(w: &mut Workflow) -> InsertSizer {
        w.cfg.set_bool_env(&[IS_COMPOSITE_GENOME]);
        w.cfg.set_string_env(&[FILTERED_INSERT_SIZES_FILE, FILTERED_STEM_LENGTHS_FILE]);

        cfg.set_bool(IS_SIZE_SELECTED, cfg.get_usize(MIN_ALLOWED_SIZE) > 0 || 
                                                cfg.get_usize(MIN_SELECTED_SIZE) > 0);
        if cfg.get_usize(MIN_ALLOWED_SIZE) == 0 {
            cfg.set_usize(MIN_ALLOWED_SIZE, (*cfg.get_usize(MIN_SELECTED_SIZE) as f64 / 
                                                        (1.0 + cfg.get_f64(SELECTED_SIZE_CV))) as usize);
        }
        cfg.set_usize(MAX_ALLOWED_SIZE, cfg.get_usize(MIN_ALLOWED_SIZE) * 2);
        cfg.set_usize(SIZE_PLOT_BIN_SIZE, if cfg.equals_string(READ_LENGTH_TYPE, "short") { 10 } else { 250 });


        w.ctrs.add_indexed_counters(&[
            (INSERT_SIZES_CHIMERIC,     "distribution of chimeric insert sizes"),
            (INSERT_SIZES_NON_CHIMERIC, "distribution of non-chimeric insert sizes"),
            (STEM_LENGTHS_CHIMERIC,     "distribution of chimeric junction stem lengths"),
            (STEM_LENGTHS_NON_CHIMERIC, "distribution of non-chimeric junction stem lengths"),
        ]);
        InsertSizer{
            is_composite_genome: *w.cfg.get_bool(IS_COMPOSITE_GENOME),
            translocation_types: vec![
                format!("{}{}{}", INTRAGENOMIC, CHIMERIC_DELIMITER, CHIMERIC),
                format!("{}{}{}", INTRAGENOMIC, CHIMERIC_DELIMITER, NOT_CHIMERIC),
                format!("{}{}{}", INTERGENOMIC, CHIMERIC_DELIMITER, CHIMERIC),
                format!("{}{}{}", INTERGENOMIC, CHIMERIC_DELIMITER, NOT_CHIMERIC),
            ],
            insert_size: 0,
            stem_length: 0,
        }
    }

    /* ---------------------------------------------------------------------------
    export likely SV artifact reads for exploring artifact mechanisms
    ---------------------------------------------------------------------------- */
    /// Record likely SV artifact reads for printing to a separate file.
    /// Only on-target reads with exactly two alignments reach this sub.
    pub fn record_variant_sizes(
        &self, 
        w: &mut Workflow, 
        tool: &mut Tool,
        aln1: SamRecord, 
        aln2: SamRecord,
        jxn_failure_flag: JxnFailureFlag,
    ) {
        // do not record foldback inversions as intermolecular chimeras
        match jxn_failure_flag {
            JxnFailureFlag::FoldbackInv => return,
            _ => {}
        }

        // let stem_length = 
    //     my $stemLengthBin = int(min(
    //         $alns[0][STEM5_LENGTH],
    //         $alns[1][STEM3_LENGTH] ? $alns[1][STEM3_LENGTH] : 1e9
    //     ) / $sizePlotBinSize) * $sizePlotBinSize;

        // record insert sizes for all single-junction SV reads stratified by chimericity
        // non-chimeric here may contain deletion and other true SVs
        let is_chimeric = is_chimeric_jxn(jxn_failure_flag);
        if is_chimeric {
            w.ctrs.increment_indexed(INSERT_SIZES_CHIMERIC, self.insert_size);
            w.ctrs.increment_indexed(STEM_LENGTHS_CHIMERIC, self.stem_length);
            // chimeric
        } else {
            w.ctrs.increment_indexed(INSERT_SIZES_NON_CHIMERIC, self.insert_size);
            w.ctrs.increment_indexed(STEM_LENGTHS_NON_CHIMERIC, self.stem_length);
            // non-chimeric
        }

        // only print single-junction interchromosomal chimeras as presumptive artifacts
        // we do not yet know if these are single-molecule or would be confirmed downstream
        // but most interchromosomal, i.e., translocation molecules are likely artifacts
        if aln1.rname == aln2.rname { return; }

        // stratify interchromosomal molecules as inter- or intra-genomic, when applicable
        let genomicity = if tool.chroms.is_same_genome_suffix(&aln1.rname, &aln2.rname) {
            INTRAGENOMIC
        } else {
            INTERGENOMIC
        };

        // reads were stratified above as chimeric or not
        // "chimeric" means a single junction that failed breakpointMatchesSite, isOntFollowOn, or junctionHasAdapters
        // "not chimeric" means it passed those chimeric tests, but it may still be chimeric by some other mechanism
        let chimericity = if is_chimeric { CHIMERIC } else { NOT_CHIMERIC };
        let chimera_type = format!("{}{}{}", genomicity, CHIMERIC_DELIMITER, chimericity);

        // increment size counts for translocation presumptive artifacts
        // $insertSizeCounts{$type}{$sizeBin}++;
        // $stemLengthCounts{$type}{$stemLengthBin}++;

    }

    /// print summary information
    pub fn print_chimera_summary(){}
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

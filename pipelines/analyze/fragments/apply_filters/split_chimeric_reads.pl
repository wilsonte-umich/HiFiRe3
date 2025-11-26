use strict;
use warnings;

# action:
#   split sequences/molecules based on node distance to closest RE site
#       junctions are rejected if either node is <= REJECT_JUNCTION_DISTANCE bp from an RE site
#       the flanking alignments are consistent with a foldback inversion 
#       adapter is present in junction inserted bases of sufficient length
#   if ONT, additionally split sequences/molecules to contigs if
#       junction inserted bases have a low quality span, e.g., as results from an ONT follow-on event
#   when a junction is rejected as chimeric, return an appropriate flag to suppress reporting of that junction downstream

#----------------------------------------------------------------------------
# initialize
#----------------------------------------------------------------------------
# initialize reporting
our (@siteDistances, @adapterScores, %insertSizeCounts, %stemLengthCounts);

# variables
use vars qw(
    $REJECTING_JUNCTION_RE_SITES
    $REJECT_JUNCTION_DISTANCE
    $INSERTION_WINDOW_SIZE
    $MIN_INSERTION_WINDOW_QUAL
    $INSERTION_ADAPTER_SEQUENCE
    $MIN_ADAPTER_LENGTH
    $MAX_ADAPTER_LENGTH
    $MIN_ADAPTER_SCORE
    $IS_COMPOSITE_GENOME
    $FILTERED_INSERT_SIZES_FILE
    $FILTERED_STEM_LENGTHS_FILE
    %chromIndex
    $isONT
    $sizePlotBinSize
);
my $minSumInsQual = ($MIN_INSERTION_WINDOW_QUAL + 33) * $INSERTION_WINDOW_SIZE;
my $adapterCore = $INSERTION_ADAPTER_SEQUENCE; # duplex portion of the ONT kit adapter; for ligation kit, last T matches the one-base A-tail; fused to 5' genomic ends 
my $adapterCoreRc = $adapterCore;    # fused to 3' genomic ends
rc(\$adapterCoreRc);
my ($aln1, $aln2, $readLen,
    $isChimeric, $isFoldback, $clip1, $clip2, $jxnInsSize);

# constants
use constant {
    S_QNAME             => 0, # SITE_SAM fields
    S_FLAG              => 1,
    S_RNAME             => 2,
    S_POS1              => 3,
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    SITE_INDEX1_1       => 6,
    SITE_POS1_1         => 7,
    SITE_DIST_1         => 8,
    SITE_INDEX1_2       => 9,
    SITE_POS1_2         => 10,
    SITE_DIST_2         => 11,
    SEQ_SITE_INDEX1_2   => 12,
    SEQ_SITE_POS1_2     => 13,
    IS_END_TO_END_READ  => 14,
    IS_END_TO_END_INSERT=> 15,
    NODE_5              => 16,
    NODE_3              => 17,
    CH_TAG              => 18,
    TL_TAG              => 19,
    INSERT_SIZE         => 20,
    IS_ALLOWED_SIZE     => 21,
    FM_TAG              => 22,
    DE_TAG              => 23,
    HV_TAG              => 24,
    N_REF_BASES         => 25,
    N_READ_BASES        => 26,
    STEM5_LENGTH        => 27,
    STEM3_LENGTH        => 28,
    # ... unused fields ...
    S_SEQ               => 37,
    S_QUAL              => 38,
    CS_TAG              => 39,
    #-------------
    _IS_PAIRED      => 1,  # SAM FLAG bits
    _PROPER_PAIR    => 2,
    _UNMAPPED       => 4,
    _MATE_UNMAPPED  => 8,
    _REVERSE        => 16,
    _MATE_REVERSE   => 32,
    _FIRST_IN_PAIR  => 64,
    _SECOND_IN_PAIR => 128,
    _SECONDARY      => 256,
    _FAILED_QC      => 512,
    _DUPLICATE      => 1024,
    _SUPPLEMENTAL   => 2048,
    #--------------
    MAX_SITE_DISTANCE => 100, # just used for tabulating junction site distances, not a filtering threshold
    #--------------
    FALSE => 0,
    TRUE  => 1,
    # -------------
    JXN_FAIL_NONE            => 0,
    JXN_FAIL_TRAVERSAL_DELTA => 1,
    JXN_FAIL_NONCANONICAL    => 2,
    JXN_FAIL_SITE_MATCH      => 4,
    JXN_FAIL_FOLDBACK_INV    => 8,
    JXN_FAIL_ONT_FOLLOW_ON   => 16, 
    JXN_FAIL_HAS_ADAPTER     => 32,
    # -------------
    READ_LEN_PROPER       => "nonSV.actual",
    PROJ_LEN_PROPER       => "nonSV.projected",
    READ_LEN_CHIMERIC     => "SV.chimeric", # cannot assess projected size of structural variant reads
    READ_LEN_NON_CHIMERIC => "SV.passed",   # these include all SV reads, not just translocations
    #-------------
    INTRAGENOMIC  => "intra", # for assembling types of translocation, i.e., presumed artifacts
    INTERGENOMIC  => "inter",
    CHIMERIC      => "chimeric",
    NOT_CHIMERIC  => "passed",
    CHIMERIC_DELIMITER => ".",
};
my @translocationTypes = (
    INTRAGENOMIC.CHIMERIC_DELIMITER.CHIMERIC, 
    INTRAGENOMIC.CHIMERIC_DELIMITER.NOT_CHIMERIC, 
    INTERGENOMIC.CHIMERIC_DELIMITER.CHIMERIC, 
    INTERGENOMIC.CHIMERIC_DELIMITER.NOT_CHIMERIC
);
my $READ_LEN_CHIMERIC = READ_LEN_CHIMERIC;
my $READ_LEN_NON_CHIMERIC = READ_LEN_NON_CHIMERIC;

# check junction quality from pairs of alignments
sub getJxnFailureFlag {
    ($aln1, $aln2, $readLen) = @_;
    $isChimeric = FALSE;
    $isFoldback = FALSE;

    # reject junction with an alignment on a non-canonical chromosome
    $chromIndex{$$aln1[S_RNAME]} or return JXN_FAIL_NONCANONICAL;
    $chromIndex{$$aln2[S_RNAME]} or return JXN_FAIL_NONCANONICAL;

    # reject junctions that fail REJECT_JUNCTION_DISTANCE
    # reject ONT junctions that contain adapters, if not already caught by REJECT_JUNCTION_DISTANCE    
    $REJECTING_JUNCTION_RE_SITES and breakpointMatchesSite() and return JXN_FAIL_SITE_MATCH;
    isFoldbackInversion()   and return JXN_FAIL_FOLDBACK_INV;
    isOntFollowOn()         and return JXN_FAIL_ONT_FOLLOW_ON;
    junctionHasAdapters()   and return JXN_FAIL_HAS_ADAPTER;

    # junction passed all criteria
    return JXN_FAIL_NONE; 
}
sub setChimeric {
    $isChimeric = TRUE;
    TRUE;
}

#----------------------------------------------------------------------------
# chimera filter functions
#----------------------------------------------------------------------------
# check whether junction breakpoints are too close to a known RE filtering site
sub breakpointMatchesSite {
    my $dist1 = abs( $$aln1[SITE_DIST_2] ); # junction node 1, at query 3' end of left  flanking alignment
    my $dist2 = abs( $$aln2[SITE_DIST_1] ); # junction node 2, at query 5' end of right flanking alignment
    $siteDistances[min($dist1, MAX_SITE_DISTANCE)]++;
    $siteDistances[min($dist2, MAX_SITE_DISTANCE)]++;
    min($dist1, $dist2) > $REJECT_JUNCTION_DISTANCE and return FALSE;
    return setChimeric();
}
# check whether a sequence is consistent with an ONT duplex read as a single foldback inversion
# it is most sensitive and acceptable to reject any inversion junction with reverse-complement overlap between its flanking alignments
# UPDATE: do this for all libraries now, not just ONT
#   ----->
#         | inversion junction
#   <-----
sub isFoldbackInversion {
    $$aln1[S_RNAME] eq $$aln2[S_RNAME] or return FALSE; # translocation
    ($$aln1[S_FLAG] & _REVERSE) != ($$aln2[S_FLAG] & _REVERSE) or return FALSE; # deletion or duplication
    (
        $$aln1[S_POS1] <= getEnd(@{$aln2}[S_POS1, S_CIGAR]) and
        $$aln2[S_POS1] <= getEnd(@{$aln1}[S_POS1, S_CIGAR])
    ) or return FALSE;
    $isFoldback = TRUE;
    return TRUE; # do not increment chimeric size count for foldback inversions, they are intra, not inter-molecular
}

# determine if an ONT junction has a very low quality insertion span, identifying it as a two-insert follow-on event
sub isOntFollowOn {

    # locate the insertion
    # stop if insertion of insufficient size to warrant adapter/follow-on detection
    # do this in advance of junctionHasAdapters even if not an ONT library
    $clip1 = 
        ($$aln1[S_FLAG] & _REVERSE) ? 
        getLeftClip($$aln1[S_CIGAR]) :
        getRightClip($$aln1[S_CIGAR]);
    $clip2 = 
        ($$aln2[S_FLAG] & _REVERSE) ? 
        getRightClip($$aln2[S_CIGAR]) :
        getLeftClip($$aln2[S_CIGAR]);
    $jxnInsSize = ($clip1 + $clip2) - $readLen;
    $isONT or return FALSE;
    $jxnInsSize < $MIN_ADAPTER_LENGTH and return FALSE;
    $jxnInsSize > $INSERTION_WINDOW_SIZE or return FALSE;

    # examine the inserted bases for low-quality base stretches as occurs during follow-on missing signals
    my $insertedQuals = 
        ($$aln1[S_FLAG] & _REVERSE) ? 
        substr(substr($$aln1[S_QUAL], 0, $clip1), -$jxnInsSize) :
        substr(substr($$aln1[S_QUAL], -$clip1), 0, $jxnInsSize);
    my @insertedQuals = split("", $insertedQuals);
    my $sumInsQual = 0;
    map{ $sumInsQual += ord($_) } @insertedQuals[0..($INSERTION_WINDOW_SIZE - 1)];
    if($sumInsQual < $minSumInsQual){
        # $nInsertionCuts++;
        return setChimeric();
    }
    # examine windows of bases to find no-base signal stretches
    # observed behavior is read1 ... (3' adapter) ... no-base signals/bases ... 5' adapter ... read2
    foreach my $i(1..($jxnInsSize - $INSERTION_WINDOW_SIZE)){
        $sumInsQual -= ord($insertedQuals[$i - 1]);
        $sumInsQual += ord($insertedQuals[$i + $INSERTION_WINDOW_SIZE - 1]);
        if($sumInsQual < $minSumInsQual){
            # $nInsertionCuts++;
            return setChimeric();
        }
    }
    return FALSE;
}

# determine if an ONT junction has adapters, identifying it as a two-insert event
sub junctionHasAdapters {
    $jxnInsSize < $MIN_ADAPTER_LENGTH and return FALSE;
    $jxnInsSize > $MAX_ADAPTER_LENGTH and return FALSE;

    # next examine the inserted bases for ONT adapters using Smith-Waterman on both strands
    my $insertedBases = 
        ($$aln1[S_FLAG] & _REVERSE) ? 
        substr(substr($$aln1[S_SEQ], 0, $clip1), -$jxnInsSize) :
        substr(substr($$aln1[S_SEQ], -$clip1), 0, $jxnInsSize);
    my ($qryOnRef1, $score1) = smith_waterman($insertedBases, $adapterCore);
    $adapterScores[max(0, min($score1, 100))]++;
    if($score1 >= $MIN_ADAPTER_SCORE){
        return setChimeric();
    }
    my ($qryOnRef2, $score2) = smith_waterman($insertedBases, $adapterCoreRc);
    $adapterScores[max(0, min($score2, 100))]++;
    if($score2 >= $MIN_ADAPTER_SCORE){
        return setChimeric();
    }
    return FALSE;
}

#----------------------------------------------------------------------------
# export likely SV artifact reads for exploring artifact mechanisms
#----------------------------------------------------------------------------
# export likely SV artifact reads to a separate file
# only on-target reads with exactly two alignments reach this sub
sub recordVariantSizes {
    my ($sizeBin, @alns) = @_;

    # do not record foldback inversions as intermolecular chimeras
    $isFoldback and return;
    my $stemLengthBin = int(min(
        $alns[0][STEM5_LENGTH],
        $alns[1][STEM3_LENGTH] ? $alns[1][STEM3_LENGTH] : 1e9
    ) / $sizePlotBinSize) * $sizePlotBinSize;

    # record insert sizes for all single-junction SV reads stratified by chimericity
    # non-chimeric here may contain deletion and other true SVs
    if($isChimeric){
        $insertSizeCounts{$READ_LEN_CHIMERIC}{$sizeBin}++;
        $stemLengthCounts{$READ_LEN_CHIMERIC}{$stemLengthBin}++;
    } else {
        $insertSizeCounts{$READ_LEN_NON_CHIMERIC}{$sizeBin}++;
        $stemLengthCounts{$READ_LEN_NON_CHIMERIC}{$stemLengthBin}++;
    }

    # only print single-junction interchromosomal chimeras as presumptive artifacts
    # we do not yet know if these are single-molecule or would be confirmed downstream
    # but most interchromosomal, i.e., translocation molecules are likely artifacts
    $alns[0][S_RNAME] eq $alns[1][S_RNAME] and return;

    # stratify interchromosomal molecules as inter- or intra-genomic, when applicable
    my $genomicity = $IS_COMPOSITE_GENOME ?
        isInterGenomic($alns[0][S_RNAME], $alns[1][S_RNAME]) :
        INTRAGENOMIC;

    # reads were stratified above as chimeric or not
    # "chimeric" means a single junction that failed breakpointMatchesSite, isOntFollowOn, or junctionHasAdapters
    # "not chimeric" means it passed those chimeric tests, but it may still be chimeric by some other mechanism
    my $chimericity = $isChimeric ? CHIMERIC : NOT_CHIMERIC;
    my $type = $genomicity.CHIMERIC_DELIMITER.$chimericity;

    # increment size counts for translocation presumptive artifacts
    $insertSizeCounts{$type}{$sizeBin}++;
    $stemLengthCounts{$type}{$stemLengthBin}++;
}
sub isInterGenomic {
    my ($chrom1, $chrom2) = @_;
    my ($genome1) = ($chrom1 =~ m/.+_(.+)/); # e.g., chr1_(hs1)
    my ($genome2) = ($chrom2 =~ m/.+_(.+)/);
    $genome1 eq $genome2 ? INTRAGENOMIC : INTERGENOMIC;
}

# print summary information
sub printChimeraSummary {

    # print and save insert size distributions
    open my $sizesH, '>', "$FILTERED_INSERT_SIZES_FILE" or throwError("could not open insert sizes file: $!");
    open my $stemsH, '>', "$FILTERED_STEM_LENGTHS_FILE" or throwError("could not open stem lengths file: $!");
    foreach my $type(
        READ_LEN_PROPER, PROJ_LEN_PROPER, READ_LEN_CHIMERIC, READ_LEN_NON_CHIMERIC, 
        @translocationTypes
    ){
        if($insertSizeCounts{$type}){
            my @insertSizes = sort { $a <=> $b } keys %{$insertSizeCounts{$type}};
            print $sizesH join("\n", map { 
                join("\t", 
                    $type,
                    $_, 
                    $insertSizeCounts{$type}{$_}
                )
            } @insertSizes), "\n";
        }
        if($stemLengthCounts{$type}){
            my @stemLengths = sort { $a <=> $b } keys %{$stemLengthCounts{$type}};
            print $stemsH join("\n", map { 
                join("\t", 
                    $type,
                    $_, 
                    $stemLengthCounts{$type}{$_}
                )
            } @stemLengths), "\n";
        }
    }
    close $sizesH;
    close $stemsH;

    # print site distance distribution to log
    if($REJECTING_JUNCTION_RE_SITES){
        print STDERR "\njunction node-to-site distances\n";
        foreach my $dist(0..MAX_SITE_DISTANCE){
            print STDERR join("\t", $dist, $siteDistances[$dist] || 0), "\n";
        }
    }

    # print adapter SW score distribution to log
    print STDERR "\nadapter SW scores\n";
    foreach my $score(0..100){
        print STDERR join("\t", $score, $adapterScores[$score] || 0), "\n";
    }
}

1;

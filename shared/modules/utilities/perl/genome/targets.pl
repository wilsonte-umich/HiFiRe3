use strict;
use warnings;

#----------------------------------------------------------
# manipulations related to genome target regions
#----------------------------------------------------------

# output variables
our ($nRegions, $sumTargetLens, $sumPaddedTargetLens) = (0, 0, 0);
our (%targetChroms, %targetRegions, %targetClassCounts, @regionCenters);

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
    IS_END_TO_END       => 14,
    CH_TAG              => 15,
    TL_TAG              => 16,
    INSERT_SIZE         => 17,
    IS_ALLOWED_SIZE     => 18,
    DE_TAG              => 19,
    HV_TAG              => 20,
    N_REF_BASES         => 21,
    N_READ_BASES        => 22,
    STEM5_LENGTH        => 23,
    STEM3_LENGTH        => 24,
    PASSED_STEM5        => 25,
    PASSED_STEM3        => 26,
    BLOCK_N             => 27,
    ALN_FAILURE_FLAG    => 28,
    JXN_FAILURE_FLAG    => 29,
    TARGET_CLASS        => 30,
    # ... unused fields ...
    #-------------
    ON_TARGET   => "T",
    NEAR_TARGET => "A", # for "adjacent"
    OFF_TARGET  => "-",
    #-------------
    NULL_TARGET_I    => 0,
    NULL_TARGET_NAME => "*",
    NULL_DISTANCE    => 0,
    #-------------
    FALSE  => 0,
    TRUE   => 1,
};
my %targetClasses = (
    "--" => 0,
    "TT" => 1,
    "TA" => 2,
    "T-" => 3,
    "AA" => 4,
    "A-" => 5
);
my @countClasses = (
    OFF_TARGET,
    ON_TARGET,
    ON_TARGET,
    ON_TARGET,
    NEAR_TARGET,
    NEAR_TARGET
);

# operating parameters
use vars qw($TARGETS_BED $REGION_PADDING $TARGET_SCALAR);
defined $REGION_PADDING or $REGION_PADDING = 0;
defined $TARGET_SCALAR  or $TARGET_SCALAR  = 1;

# load the target regions
sub loadTargetRegions {
    my ($quiet) = @_;
    ($TARGETS_BED and $TARGETS_BED ne "null" and $TARGETS_BED ne "NA") or return;

    # first pass to record target regions (T=target type)
    @regionCenters = ();
    $nRegions = loadTargetRegions_(\$sumTargetLens, 0, ON_TARGET);

    # second pass to record padded regions (A=adjacent type)
    if($REGION_PADDING){
        loadTargetRegions_(\$sumPaddedTargetLens, $REGION_PADDING, NEAR_TARGET);
    } else {
        $sumPaddedTargetLens = $sumTargetLens;
    }

    # report target summary to log
    if(!$quiet){
        printCount($nRegions,            'nRegions',            'target regions');
        printCount(commify($sumTargetLens),       'sumTargetLens',       'total bp covered by target regions');
        $REGION_PADDING and 
        printCount(commify($sumPaddedTargetLens), 'sumPaddedTargetLens', 'total bp covered by padded target regions');
    }
}
sub loadTargetRegions_ {
    my ($sumLens, $regPad, $type) = @_;
    open my $inH, "<", $TARGETS_BED or die "could not open $TARGETS_BED: $!\n";
    my $targetI1 = 0;
    while(my $line = <$inH>){
        $targetI1++; # 1-referenced
        chomp $line;
        $line =~ s/\r//g;
        my ($chr, $start, $end, $name) = split("\t", $line);
        $targetChroms{$chr} = 1; # lookup for whether a chromosome has target regions
        for(my $pos =  int(($start - $regPad) / $TARGET_SCALAR);
               $pos <= int(($end   + $regPad) / $TARGET_SCALAR);
               $pos++){    
            $targetRegions{$chr} and $targetRegions{$chr}{$pos} and next;
            $targetRegions{$chr}{$pos} = [$targetI1, $type, $name]; # lookup for the target region associated with a coordinate
        }
        $$sumLens += ($end + $regPad) - ($start - $regPad);
        $type eq ON_TARGET and push @regionCenters, int(($start + $end) / 2);
    }
    close $inH;
    return $targetI1;
}
sub getTargetRegions {
    ($TARGETS_BED and $TARGETS_BED ne "null" and $TARGETS_BED ne "NA") or return [];
    my @regions;
    open my $inH, "<", $TARGETS_BED or die "could not open $TARGETS_BED: $!\n";
    while(my $line = <$inH>){
        chomp $line;
        $line =~ s/\r//g;
        my ($chr, $start, $end, $name) = split("\t", $line);
        push @regions, {
            name  => $name,
            chr   => $chr,
            start => $start,
            end   => $end,
            paddedStart  => $start - $REGION_PADDING, # all still half-open like the source BED
            paddedEnd    => $end   + $REGION_PADDING,
            paddedStart1 => $start - $REGION_PADDING + 1 # 1-referenced start, unlike above
        };
    }
    close $inH;
    return \@regions;
}
sub getTargetChroms {
    return keys %targetChroms;
}

# return information on a single genome position relative to target regions
sub getPosTarget {
    my ($chr, $pos) = @_;
    my $ct = $targetRegions{$chr} or return (NULL_TARGET_NAME, NULL_DISTANCE);
    my $ctp = $$ct{int($pos / $TARGET_SCALAR)} || 0;
    if($ctp){
        my ($targetI1, $type, $name) = @$ctp;
        ($name, $pos - $regionCenters[$targetI1 - 1]);
    } else {
        (NULL_TARGET_NAME, NULL_DISTANCE);
    }
}

# codified form of the relationship of two coordinate pairs with respect to target regions
sub getAlnTarget { 
    my ($aln) = @_;
    my $leftPos1  = $$aln[S_POS1];
    my $rightPos1 = getEnd($$aln[S_POS1], $$aln[S_CIGAR]);
    my $ct = $targetRegions{$$aln[S_RNAME]} or return (NULL_TARGET_I, NULL_TARGET_I, $leftPos1, $rightPos1); # --
    my $ct1 = $$ct{int($leftPos1  / $TARGET_SCALAR)} || 0; # check the two ends of the alignment (not the whole read)
    my $ct2 = $$ct{int($rightPos1 / $TARGET_SCALAR)} || 0;
    if($ct1 and $ct2){
        my ($targetI1_1, $type1) = @$ct1;
        my ($targetI1_2, $type2) = @$ct2;
        my $tc = join("", sort {$b cmp $a} $type1, $type2); # TT, TA or AA
        ($targetClasses{$tc}, $targetI1_1 == $targetI1_2 ? $targetI1_1 : NULL_TARGET_I, $leftPos1, $rightPos1);
    } elsif($ct1 or $ct2){
        my ($targetI1, $type) = $ct1 ? @$ct1 : @$ct2;
        ($targetClasses{$type.OFF_TARGET}, $targetI1, $leftPos1, $rightPos1); # T- or A-
    } else {
        (NULL_TARGET_I, NULL_TARGET_I, $leftPos1, $rightPos1);
    }
}

# set the target class metadata for a list of alignments
# assembled TARGET_CLASS reflects:
#   the read outermost alignments
#   each alignment, with its matching targetI1
sub setAlnTargetClasses {
    my (@alns) = @_;
    my ($targetClass5, $targetI1_5, $leftPos1_5, $rightPos1_5) = getAlnTarget($alns[0]);
    if(@alns > 1){
        my ($targetClass3) = getAlnTarget($alns[$#alns]);
        foreach my $aln(@alns){
            my ($targetClass, $targetI1) = getAlnTarget($aln);
            $$aln[TARGET_CLASS] = $targetI1 << 9 | $targetClass5 << 6 | $targetClass3 << 3 | $targetClass; 
            # only non-SV reads are used for enrichment/coverage assessment
        }
    } else {
        $alns[0][TARGET_CLASS] = $targetI1_5 << 9 | $targetClass5 << 6 | $targetClass5 << 3 | $targetClass5;
        my $cc = $countClasses[$targetClass5];
        $targetClassCounts{$cc}{nAlns}++;
        $targetClassCounts{$cc}{nSeqBases}  += $rightPos1_5 - $leftPos1_5 + 1;
        $targetClassCounts{$cc}{nProjBases} += abs($alns[0][SEQ_SITE_POS1_2] - $alns[0][SITE_POS1_1]) + 1;
    }
    return $targetI1_5 == NULL_TARGET_I ? FALSE : TRUE; # return $isOnTarget_5
}

# print target class counts to log
sub printTargetCounts_ {
    my ($cc, $N, $t, $T) = @_;
    my $tcc = $targetClassCounts{$cc} or return;
    my ($n, $r);

    $n = $$tcc{nAlns} || 0;
    $r = $n / $N;
    printCount(commify($n),          join("_", "nAlns", $t),      join(" ", "number of", $T, "alignments"));
    printCount(roundCount($r, 1000), join("_", "alnDensity", $t), join(" ", $T, "alignment density (alignments per base)"));

    $n = $$tcc{nSeqBases} || 0;
    $r = $n / $N;
    printCount(commify($n),         join("_", "nSeqBases", $t),       join(" ", "number of", $T, "sequenced bases"));
    printCount(roundCount($r, 100), join("_", "seqBaseCoverage", $t), join(" ", $T, "sequenced base fold coverage"));

    $n = $$tcc{nProjBases} || 0;
    $r = $n / $N;
    printCount(commify($n),         join("_", "nProjBases", $t),       join(" ", "number of", $T, "insert bases, including projections"));
    printCount(roundCount($r, 100), join("_", "projBaseCoverage", $t), join(" ", $T, "insert base fold coverage"));

    $r = ($$tcc{nSeqBases} || 0)  / ($$tcc{nAlns} || 0);
    printCount(roundCount($r, 100), join("_", "avgAlnLen", $t), join(" ", $T, "average alignment length"));

    $r = ($$tcc{nProjBases} || 0) / ($$tcc{nAlns} || 0);
    printCount(roundCount($r, 100), join("_", "avgInsertLen", $t), join(" ", $T, "average insert length"));

    $r = ($$tcc{nSeqBases} || 0) / ($$tcc{nProjBases} || 0);
    printCount(roundCount($r, 1000), join("_", "fracSequencedBases", $t), join(" ", $T, "fraction of insert bases sequenced"));
}
sub printTargetCounts {
    my ($GENOME_SIZE) = @_;
    my $sumPadding = $sumPaddedTargetLens - $sumTargetLens;
    my $sumOffTarget = $GENOME_SIZE - $sumPaddedTargetLens;
    printTargetCounts_(ON_TARGET,   $sumTargetLens, "ON(T)",   "on-target");
    $REGION_PADDING and 
    printTargetCounts_(NEAR_TARGET, $sumPadding,    "NEAR(A)", "near-target");
    printTargetCounts_(OFF_TARGET,  $sumOffTarget,  "OFF(-)",  "off-target");
}

1;

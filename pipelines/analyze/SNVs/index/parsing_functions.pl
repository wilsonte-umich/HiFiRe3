use strict;
use warnings;

# SNV parsing functions shared by RE and non-RE fragment parsers

# working variables
use vars qw($SNV_ALNS_PREFIX $paddedChromIndex1 $leftClip_ $rightClip_ $MIN_SNV_INDEL_QUAL);
my (@vaVarQualities);

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
    END5_ON_TARGET      => 31,
    READ_HAS_JXN        => 32,
    S_SEQ               => 33,
    S_QUAL              => 34,
    CS_TAG              => 35,
    #-------------
    S_VA_TAG            => 33,
    S_LEFT_CLIP         => 34,
    S_RIGHT_CLIP        => 35,
    S_END_POS1          => 36,
    #-------------
    FALSE  => 0,
    TRUE   => 1,
    #-------------
    CS_MATCH        => ":",
    CS_SUBSTITUTION => "*",
    CS_INSERTION    => "+",
    CS_DELETION     => "-",
    #-------------
    NULL_ENTRY => "*",
};

# open/close the required chromosome-level output files
our ($snvAlnsH);
sub openFileHandles {
    open $snvAlnsH, "|-", "gzip -c > $SNV_ALNS_PREFIX.$paddedChromIndex1.txt.gz" or die "could not open snvAlnsH: $!\n";
}
sub closeFileHandles {
    close $snvAlnsH;
}

# reverse various values for alignments matching the bottom reference strand
# this modifies the fields used downstream so that reverse alignments are now aligned 
#    to the bottom reference _strand_ as the reference coordinate system
# however, the base _values_ and variant stretches are not complemented, they still match the top reference strand
# do not flip the reverse bit in the FLAG field; it is retained as history that the alignment was flipped
sub invertReverseStrandAlignment {
    my ($aln, $chromSize) = @_;
    @$aln[S_POS1, S_END_POS1] = (
        $chromSize - $$aln[S_END_POS1] + 1,
        $chromSize - $$aln[S_POS1]     + 1
    );
    @$aln[S_LEFT_CLIP, S_RIGHT_CLIP] = @$aln[S_RIGHT_CLIP, S_LEFT_CLIP];
    $$aln[S_VA_TAG] = reverse_va_tag($$aln[S_VA_TAG]);
    $$aln[S_QUAL] = reverse($$aln[S_QUAL]);
}

# extract base qualities for variant stretches and check them against MIN_SNV_INDEL_QUAL
my $vaMatchOp_start = qr/^(\d+)(:)(\d+);/;
my $vaVarOp_start   = qr/^(\d+)(\*)([acgtn]*),([acgtn]*);/;
sub getQsTag{
    my ($aln) = @_;

    # when (haplotype-adjusted) READ_HAS_SNV is FALSE, infer that all operations in the alignment have high quality
    my $forceHighQuality = $$aln[S_QUAL] eq NULL_ENTRY ? TRUE : FALSE;

    # match va tag operations with a boolean passed-quality flag
    my $readOffset0 = $$aln[S_LEFT_CLIP]; # working left to right along reference strand
    my $vaTag = $$aln[S_VA_TAG];
    my $qsTag = "";
    while(
        $vaTag =~ s/$vaMatchOp_start// or
        $vaTag =~ s/$vaVarOp_start//
    ){
        my ($refPos1, $op, $opVal1, $opVal2) = ($1, $2, $3, $4);
        if($op eq CS_MATCH){
            $readOffset0 += $opVal1;
            $qsTag .= TRUE; # a reference match implies high base quality
        } else {
            my $nAltBases = length($opVal2);
            my $avgBaseQual;
            # if($opVal2 eq "n"){ # this block validates that QUAL values are correctly extracted and calculated on both strands, varying positions, etc.
            #     my $Q = substr($$aln[S_QUAL], $readOffset0, $nAltBases);
            #     print join("\t", $refPos1, $$aln[S_FLAG], $readOffset0, $op, $opVal1, $opVal2, $Q, getAvgQual($Q)), "\n";
            # }
            # 43668766        16      73      *       a       n       !       0  # n bases always have ! QUAL, so we know these are the correct positions
            # 27784876        16      73      *       a       n       !       0
            # 36579068        16      272     *       t       n       !       0
            # 47769748        0       122     *       a       n       !       0
            # 50567697        0       203     *       g       n       !       0
            # 38872950        0       34      *       a       n       !       0
            if($nAltBases > 0){ # all variants except isolated deletions, get local quality from the variant stretch itself (all bases contribute to the quality)
                $forceHighQuality or $avgBaseQual = getAvgQual(substr($$aln[S_QUAL], $readOffset0, $nAltBases));
                $readOffset0 += $nAltBases; # must do after reading the current QUAL
            } else { # isolated deletions, no alignment bases in variant, get local quality from flanking bases
                $forceHighQuality or $avgBaseQual = getAvgQual(substr($$aln[S_QUAL], $readOffset0 - 1, 2));
            }
            if($forceHighQuality){
                $qsTag .= TRUE;
            } else {
                $qsTag .= $avgBaseQual >= $MIN_SNV_INDEL_QUAL ? TRUE : FALSE;
                $vaVarQualities[int($avgBaseQual)]++;
            }
        }
    }
    $qsTag;
}

# report quality metrics
sub reportFinalMetadata_all_thread_chroms {
    my ($childN) = @_;
    open my $baseQualH, ">", "$SNV_ALNS_PREFIX.base_qualities.$childN.txt" or die "could not open baseQualH: $!\n";
    foreach my $avgBaseQual(0..$#vaVarQualities){
        print $baseQualH join("\t", $avgBaseQual, $vaVarQualities[$avgBaseQual] || 0), "\n";
    }
    close $baseQualH;
}

1;

use strict;
use warnings;

# action:
#   check whether alignments pass alignment-level quality filters:
#       all alignments, even single, non-SV alignments per read:
#           minimum MAPQ
#           maximium gap-corrected base divergence
#       alignment in SV reads only:
#           minimum alignment length, measured as reference base span
#           (ONT only) minimum average base quality across alignment
#               with Dorado, individual base qualities range from Q0 to Q50
#               see script output for the range of average base qualities observed over SV flanking alignments
#   return a bit-encoded alignment failure flag to suppress untrusted SV junctions that would include a failed alignment

# variables
our (@nAlnsRejected, @avgBaseQual);
use vars qw(
    $MIN_MAPQ
    $MAX_DIVERGENCE
    $MIN_FLANK_LEN
    $MIN_AVG_BASE_QUAL
    $HAS_BASE_ACCURACY
);

# constants
use constant {
    S_QNAME             => 0, # SITE_SAM fields
    S_FLAG              => 1,
    S_RNAME             => 2,
    S_POS1              => 3,
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    # ... unused fields ...
    S_SEQ               => 33,
    S_QUAL              => 34,
    CS_TAG              => 35,
    # -------------
    ALN_FAIL_NONE            => 0,
    ALN_FAIL_MAPQ            => 1,
    ALN_FAIL_DIVERGENCE      => 2,
    ALN_FAIL_FLANK_LEN       => 4,
    ALN_FAIL_AVG_BASE_QUAL   => 8,
    #-------------
    AVG_BASE_QUAL_BIN_SIZE => 5,
};

# check individual alignment quality
sub getAlnFailureFlag {
    my ($aln, $de, $readHasJxn) = @_;

    # rejection criteria are enforced sequentially in order of efficiency
    # i.e., frequent, easy rejections are checked first
    # later criteria are not check if an earlier criterion already failed

    # criteria enforced on all alignments, even single, non-SV alignments
    if($$aln[S_MAPQ] < $MIN_MAPQ){
        $nAlnsRejected[ALN_FAIL_MAPQ]++;
        return ALN_FAIL_MAPQ;
    }
    if($de > $MAX_DIVERGENCE){
        $nAlnsRejected[ALN_FAIL_DIVERGENCE]++;
        return ALN_FAIL_DIVERGENCE;
    }

    # criteria only enforced when reads have SV junctions, i.e., multiple alignments
    if($readHasJxn){
        if(getEnd(@$aln[S_POS1, S_CIGAR]) - $$aln[S_POS1] < $MIN_FLANK_LEN){
            $nAlnsRejected[ALN_FAIL_FLANK_LEN]++;
            return ALN_FAIL_FLANK_LEN;
        }
        if(
            !$HAS_BASE_ACCURACY and # thus, this slow check only performed on ONT or other low accuracy platform
            $$aln[S_QUAL] ne "*"
        ){
            my $avgBaseQual = getAvgQual($$aln[S_QUAL]);
            $avgBaseQual[int($avgBaseQual / AVG_BASE_QUAL_BIN_SIZE + 0.5)]++;
            if($avgBaseQual < $MIN_AVG_BASE_QUAL){
                $nAlnsRejected[ALN_FAIL_AVG_BASE_QUAL]++;
                return ALN_FAIL_AVG_BASE_QUAL;
            }
        }
    }

    # alignment passed all criteria
    $nAlnsRejected[ALN_FAIL_NONE]++;
    return ALN_FAIL_NONE; 
}

1;

use strict;
use warnings;

# action:
#   convert a cs tag (the minimap2 _difference_ tag) to a va tag (a _variant_ tag)
#   a cs tag is a string of operations that describe the differences between an alignment and a reference
#   the cs format is defined by minimap2 as:
#       /(:\d+|\*[acgtn][acgtn]|[\+\-][acgtn]+)/g, where:
#           :\d+              identical, i.e., matching, sequence length
#           \*[acgtn][acgtn]  substitution of ref to query, one substituted base at a time
#           \+[acgtn]+        insertion to the reference
#           \-[acgtn]+        deletion from the reference
#   a va tag is a reformatting of cs operations that describes contiguous variant regions of an alignment
#       where all sequential cs variant operations are reduced to a single va * substitution operation
#   the va tag is used to create a pileup where a base can declare:
#       va: "I am part of this specific variant stretch", not simply
#       cs: "I differ from the reference at this base in this way", which doesn't relate adjacent differences to each other
#   a va tag is encoded as:
#       /\d+(:\d+|\*[acgtn]*,[acgtn]*);/g
#   where in each concatenated operation unit:
#       \d+ is the leftmost reference coordinate of the matching or variant stretch
#       :\d+ is the same as the cs tag, i.e, an identical sequence length
#       \*[acgtn]*,[acgtn]* is a single substitution operation that defines a contiguous variant stretch, where:
#           the first  [acgtn]* is zero to many reference bases, and
#           the second [acgtn]* is zero to many read alternative bases that replace the reference bases
#   properties of the va tag include:
#       all operation blocks end with a semicolon, even the last one
#       * operations are intuitively read as "these contiguous reference bases were replaced by these contiguous alternative bases"
#           substitution operations act as a UID for a contiguous SNV or indel variant
#           va tags act as a UID for a (portion of a) variant allele
#   as used in HiFiRe3, the following are additionally true:
#       terminal * + - cs operations are removed before constructing the va tag as they may not describe complete variants, thus:
#           HiFiRe3 va tags always begins and end with : operations, and
#           a va * operation is always flanked by anchoring : operations
#       va * operations define a complete variant stretch as observed in one read
#           note that sequencing errors may cause the same true variant to have different signatures in different reads 

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
    S_VA_TAG            => 27,
    S_LEFT_CLIP         => 28,
    S_RIGHT_CLIP        => 29,
    S_END_POS1          => 30,
};

# convert a cs tag to a va tag
my $csOp_start = qr/^(:\d+|\*[acgtn][acgtn]|[\+\-][acgtn]+)/;
my $csMatch_start = qr/^:(\d+)/;
my $csSub_start   = qr/^\*(.)(.)/;
my $csIns_start   = qr/^\+([acgtn]+)/;
my $csDel_start   = qr/^\-([acgtn]+)/;
sub cs_to_va_tag {
    my ($csTag, $refPos1) = @_;
    my $refBases = "";
    my $altBases = "";  
    my $vaTag = "";

    # work one cs operation at a time
    while($csTag =~ s/$csOp_start//){
        my $csOpVal = $1;

        # for match : operations...
        if($csOpVal =~ /$csMatch_start/){

            # ...first commit any prior variant stretch
            if($refBases or $altBases){
                $vaTag .= $refPos1."*".$refBases.",".$altBases.";";
                $refPos1 += length($refBases);
                $refBases = $altBases = "";
            }

            # ...then commit the match : operation
            $vaTag .= $refPos1.$csOpVal.";";
            $refPos1 += $1;
        
        # for variant * + - operations, concatenate the reference and alignment base spans
        } elsif($csOpVal =~ /$csSub_start/){
            $refBases .= $1;
            $altBases .= $2;
        } elsif($csOpVal =~ /$csIns_start/){
            $altBases .= $1;
        } elsif($csOpVal =~ /$csDel_start/){
            $refBases .= $1;
        }
    }

    # return the final va tag
    $vaTag;
}

# convert a cs tag to a va tag in the context of analyze/SNVs indexing
# additionally sets S_LEFT_CLIP, S_RIGHT_CLIP, and S_END_POS1
# trims terminal * + - operations, and adjust all clips and positions accordingly
# this function is called before alignment reversal on the bottom strand
use vars qw($leftClip_ $rightClip_);
my $csVarOp_start = qr/^([\*\+\-])([acgtn]+)/;
my $csVarOp_end   =  qr/([\*\+\-])([acgtn]+)$/;
sub cs_to_va_tag_aln {
    my ($aln, $alnHasClipsSet) = @_;

    # set initial values of S_LEFT_CLIP, S_RIGHT_CLIP, and S_END_POS1
    # S_CIGAR is not used for indexing after this code block executes
    # does not apply to PAF-formatted alignments provided by `genotype SNVs`, which presets these values
    unless($alnHasClipsSet){
        $$aln[S_LEFT_CLIP]  = $$aln[S_CIGAR] =~ m/$leftClip_/  ? $1 : 0;
        $$aln[S_RIGHT_CLIP] = $$aln[S_CIGAR] =~ m/$rightClip_/ ? $1 : 0;
        $$aln[S_END_POS1]   = getEnd(@$aln[S_POS1, S_CIGAR]);
    }

    # remove terminal * + - operations before constructing the va tag
    # adjust clips and positions to account for the removed bases, i.e., to narrow the called alignment
    # so that variant stretches in the final va tag are always flanked by anchoring reference matches
    my $csTag = $$aln[CS_TAG];
    while($csTag =~ s/$csVarOp_start//){
        if($1 eq "*"){
            $$aln[S_LEFT_CLIP]++;
            $$aln[S_POS1]++;
        } elsif($1 eq "+"){
            $$aln[S_LEFT_CLIP] += length($2);
        } elsif($1 eq "-"){
            $$aln[S_POS1] += length($2);
        }
    }
    while($csTag =~ s/$csVarOp_end//){
        if($1 eq "*"){
            $$aln[S_RIGHT_CLIP]++;
            $$aln[S_END_POS1]--;
        } elsif($1 eq "+"){
            $$aln[S_RIGHT_CLIP] += length($2);
        } elsif($1 eq "-"){
            $$aln[S_END_POS1] -= length($2);
        }
    }

    # convert the trimmed cs tag to a va tag
    $$aln[S_VA_TAG] = cs_to_va_tag($csTag, $$aln[S_POS1]);
}

# reverse the order of operations in a va tag to support the HiFiRe3 pileup strategy
# the content of each individual va tag operation is not reversed only their order on the alignment
# thus, variant stretches on the top and bottom strands retain the same variant type UID
my $vaMatchOp_start = qr/^(\d+:\d+;)/;
my $vaVarOp_start   = qr/^(\d+\*[acgtn]*,[acgtn]*;)/;
sub reverse_va_tag {
    my ($vaTag) = @_;
    my @operations;
    while(
        $vaTag =~ s/$vaMatchOp_start// or
        $vaTag =~ s/$vaVarOp_start//
    ){
        unshift @operations, $1;
    }
    join('', @operations);
}

1;

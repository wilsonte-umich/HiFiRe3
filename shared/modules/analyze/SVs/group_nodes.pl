use strict;
use warnings;

# environment variables
fillEnvVar(\our $GENOME_FASTA,              'GENOME_FASTA');
fillEnvVar(\our $GROUP_BREAKPOINT_DISTANCE, 'GROUP_BREAKPOINT_DISTANCE');
fillEnvVar(\our $GROUP_STEM_DISTANCE,       'GROUP_STEM_DISTANCE');

# constants
use constant {
    OJXN_CHROM_INDEX1_1     => 0, # input columns
    OJXN_REF_POS1_1         => 1,
    OJXN_STRAND_INDEX0_1    => 2,
    OJXN_CHROM_INDEX1_2     => 3,
    OJXN_REF_POS1_2         => 4,
    OJXN_STRAND_INDEX0_2    => 5,
    OJXN_JXN_TYPE           => 6,
    UJXN_N_OBSERVED         => 7,
    UJXN_ALN_OFFSET         => 8,
    # ...unused fields...
    #-----------------
    DISCARD_JUNCTION => -1,
    #-----------------
    LEFTWARD      => 0,
    RIGHTWARD     => 1,
    TOP_STRAND    => 0,
    BOTTOM_STRAND => 1,
};

# variables
use vars qw($nInputJxns $nOutputJxns);

# initialize the genome and output files
setCanonicalChroms();
our @chromSizes = getChromIndexSizes("$GENOME_FASTA.fai");

# assess whether to start a new node1 group
our ($prevLocusKey_1, $prevRefPos1_1, @jxns_1, @jxns_2, @jxns_3);
sub findJunctionGroups {
    while (my $jxn = <STDIN>){
        $nInputJxns++;
        chomp $jxn;
        my @jxn = split("\t", $jxn);
        my $locusKey_1 = join(":", @jxn[OJXN_CHROM_INDEX1_1, OJXN_STRAND_INDEX0_1]);
        if($prevLocusKey_1 and 
        (
            $prevLocusKey_1 ne $locusKey_1 or 
            $prevRefPos1_1 < $jxn[OJXN_REF_POS1_1] - $GROUP_BREAKPOINT_DISTANCE # establish node1 position breaks
        )
        ){
            processJunctions_1();
            @jxns_1 = @jxns_2 = @jxns_3 = ();
        }
        push @jxns_1, \@jxn;
        $prevLocusKey_1 = $locusKey_1;
        $prevRefPos1_1 = $jxn[OJXN_REF_POS1_1];
    }
    processJunctions_1();
}

# reparse a set of junctions that fuzzy-matched at node1 to establish node2 subgroups
sub processJunctions_1 {
    if(@jxns_1 > 1){
        my ($prevLocusKey_2, $prevRefPos1_2);
        @jxns_2 = @jxns_3 = ();
        foreach my $jxn(sort {
            $$a[OJXN_CHROM_INDEX1_2]  <=> $$b[OJXN_CHROM_INDEX1_2]  or 
            $$a[OJXN_STRAND_INDEX0_2] <=> $$b[OJXN_STRAND_INDEX0_2] or
            $$a[OJXN_REF_POS1_2]      <=> $$b[OJXN_REF_POS1_2]
        } @jxns_1){
            my $locusKey_2 = join(":", @$jxn[OJXN_CHROM_INDEX1_2, OJXN_STRAND_INDEX0_2]);
            if($prevLocusKey_2 and 
              (
                $prevLocusKey_2 ne $locusKey_2 or 
                $prevRefPos1_2 < $$jxn[OJXN_REF_POS1_2] - $GROUP_BREAKPOINT_DISTANCE # establish node2 position breaks
              )
            ){
                processJunctions_2();
                @jxns_2 = @jxns_3 = ();
            }
            push @jxns_2, $jxn;
            $prevLocusKey_2 = $locusKey_2;
            $prevRefPos1_2 = $$jxn[OJXN_REF_POS1_2];
        }
        processJunctions_2();
    } else {
        print join("\t", @{$jxns_1[0]}), "\n";
        $nOutputJxns++;
    }
}

# reparse a set of junctions that fuzzy-matched at node1 and node2 to establish adjStemSum subgroups
sub processJunctions_2 {
    if(@jxns_2 > 1){
        my @adjStemSums = map { getAdjStemSum($_) } @jxns_2;
        @jxns_3 = ();
        foreach my $jxn2I(sort {
            $adjStemSums[$a] <=> $adjStemSums[$b]
        } 0..$#jxns_2){
            if(
                $jxn2I > 0 and 
                $adjStemSums[$jxn2I - 1] < $adjStemSums[$jxn2I] - $GROUP_STEM_DISTANCE # establish breaks
            ){
                processJunctions_3();
                @jxns_3 = ();
            }
            push @jxns_3, $jxns_2[$jxn2I];
        }
        processJunctions_3();
    }
    foreach my $jxn(@jxns_2){
        $$jxn[0] == DISCARD_JUNCTION and next; 
        print join("\t", @$jxn), "\n";
        $nOutputJxns++;
    }
}

# estimate the number of genome bp still present as the chromosome stem length on each side of the junction
# 1     7                   9       1 (i.e., count rightward stems from the end of the chrom)
# -----------------------------------
#    1--> (7 bp stem)       <--1 (9 bp stem)
#    <--2 (7 bp stem)       2--> (9 bp stem)
# inappropriately clipped junction bases always do two things to ~the same degree:
#   decrease the summed stem length
#   increase the size of the alignment offset, i.e., the number of apparently inserted bases
# thus, summed stem length + alignment offset = adjusted stem sum stays ~constant as a clip-resistant metric
# adjusted stem sum is a junction/SV-level property, unlike breakpoint distances which are assessed per breakpoint
# for efficiency, each breakpoint is assessed more promiscuously in series to establish candidate junction groups
# then the group is subjected to stem sum matching more stringently since clip errors can now be accounted for
# it is not practical to account for clip errors during breakpoint distance matching
# because we cannot know which breakpoint (either, both, neither) may have been inappropriately clipped
# adjusted stem sum, in considering both breakpoints at once, doesn't care which breakpoint was inappropriately clipped
# the resulting value is the same for:
#    no inappropriate clipping
#    read 1 inappropriately clipped
#    read 2 inappropriately clipped
#    both read 1 and read 2 inappropriately clipped
# demonstrating the resilience of this metric
# one cannot only use adjusted stem sums since wildly different junction positions could yield the same value
# but once breakpoints are known to be roughly co-localized, adjusted stem sums provide robust junction matching
# the method is robust to inversions and all types of translocations
sub getAdjStemSum {
    my $jxn = shift;
    getStemLength(
        $$jxn[OJXN_CHROM_INDEX1_1], 
        $$jxn[OJXN_REF_POS1_1], 
        $$jxn[OJXN_STRAND_INDEX0_1] == TOP_STRAND ? LEFTWARD : RIGHTWARD
    ) + getStemLength( # the sum of stem lengths
        $$jxn[OJXN_CHROM_INDEX1_2], 
        $$jxn[OJXN_REF_POS1_2],
        $$jxn[OJXN_STRAND_INDEX0_2] == TOP_STRAND ? RIGHTWARD : LEFTWARD
    ) + $$jxn[UJXN_ALN_OFFSET] # the adjustment for alignment offset
} 
sub getStemLength {
    my ($chromIndex1, $refPos1, $side0) = @_;
    $side0 == LEFTWARD ? 
        $refPos1 : 
        $chromSizes[$chromIndex1] - $refPos1 + 1
}

1;

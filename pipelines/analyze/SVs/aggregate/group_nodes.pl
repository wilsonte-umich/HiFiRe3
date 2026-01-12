use strict;
use warnings;

# actions:
#   perform fuzzy matching of junction nodes to each other to aggregate inexact junction matches
# input:
#   unique junctions with merged source metadata on STDIN
# output:
#   unique junctions now further aggregated by fuzzy node matching on STDOUT

# initialize reporting
our $script = "group_nodes";
our $error  = "$script error";
our ($nInputJxns, $nOutputJxns) = (0) x 10;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$ENV{MODULES_DIR}/analyze/SVs/$_.pl" } qw(group_nodes);
resetCountFile();

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
    UJXN_JXN_BASES          => 9,
    UJXN_PATHS              => 10,
    UJXN_N_PATH_JUNCTIONS   => 11,
    UJXN_JSRC_MAPQ          => 12,
    UJXN_JSRC_DE_TAG        => 13,
    UJXN_JSRC_SITE_DIST     => 14,
    UJXN_JSRC_ALN_FAILURE_FLAG => 15,
    UJXN_JSRC_JXN_FAILURE_FLAG => 16,
    UJXN_JSRC_QNAMES           => 17,
    UJXN_JSRC_INSERT_SIZES     => 18,
    UJXN_JSRC_IS_ALLOWED_SIZES => 19,
    UJXN_JSRC_MIN_STEM_LENGTHS => 20,
    UJXN_JSRC_PASSED_STEMS     => 21,
    UJXN_JSRC_SEQS            => 22,
    UJXN_JSRC_QUALS           => 23,
    UJXN_JSRC_CIGARS          => 24,
    UJXN_JSRC_ORIENTATIONS    => 25,
    UJXN_JSRC_CHANNELS        => 26, # concatenated lists that are dropped later after used for duplicate purging
    UJXN_JSRC_OUTER_ENDPOINTS => 27,
    #-----------------
    DISCARD_JUNCTION => -1,
};

# variables
use vars qw(@jxns_3);

# process junctions one key at a time
findJunctionGroups();

# print summary information
printCount(commify($nInputJxns),  'nInputJxns',  'input junctions');
printCount(commify($nOutputJxns), 'nOutputJxns', 'output junctions after fuzzy node matching');

# collapse a set of junctions that fuzzy-matched at node1, node2, and adjStemSum to a single output junction with adjusted metadata
sub processJunctions_3 {

    # no fuzzy matches to aggregate, a single-junction group
    @jxns_3 > 1 or return; 

    # keep the mostly highly represented junction in the subgroup as presumptive true nodes
    # node1,node2,jxnType,alnOffset,jxnBases all persist as is for this best junction (values from other jxns dropped)
    @jxns_3 = sort { $$b[UJXN_N_OBSERVED] <=> $$a[UJXN_N_OBSERVED] } @jxns_3; 

    # UJXN_N_OBSERVED summed over all aggregating junctions
    foreach my $discardedJxn(@jxns_3[1..$#jxns_3]){
        $jxns_3[0][UJXN_N_OBSERVED] += $$discardedJxn[UJXN_N_OBSERVED];

        # concatenate most values onto the best junction's aggregated data to retain all original read data
        # all fields made comma-safe upstream, e.g., MM and ML tags use : (colon) as a separator, etc.
        foreach my $i(
            UJXN_PATHS, 
            UJXN_JSRC_QNAMES, 
            UJXN_JSRC_INSERT_SIZES, UJXN_JSRC_IS_ALLOWED_SIZES, UJXN_JSRC_MIN_STEM_LENGTHS, UJXN_JSRC_PASSED_STEMS,
            UJXN_JSRC_SEQS, UJXN_JSRC_QUALS, UJXN_JSRC_CIGARS,
            UJXN_JSRC_ORIENTATIONS, UJXN_JSRC_CHANNELS, UJXN_JSRC_OUTER_ENDPOINTS,
            UJXN_JSRC_ALN_FAILURE_FLAG, UJXN_JSRC_JXN_FAILURE_FLAG # or move to best-worst min below?
        ){
            $jxns_3[0][$i] .= ",$$discardedJxn[$i]";
        }

        # aggregate UJXN_N_PATH_JUNCTIONS, UJXN_JSRC_MAPQ, UJXN_JSRC_SITE_DIST, UJXN_READ_HAS_SV, UJXN_JSRC_DE_TAG 
        #   to the "best-worst" value over all group junctions
        foreach my $i(UJXN_N_PATH_JUNCTIONS, UJXN_JSRC_MAPQ, UJXN_JSRC_SITE_DIST){
            $jxns_3[0][$i] = max($jxns_3[0][$i], $$discardedJxn[$i]);
        }
        foreach my $i(UJXN_JSRC_DE_TAG){
            $jxns_3[0][$i] = min($jxns_3[0][$i], $$discardedJxn[$i]);
        }

        # flag discarded junctions to prevent them from printing to STDOUT
        $$discardedJxn[0] = DISCARD_JUNCTION;
    }
}

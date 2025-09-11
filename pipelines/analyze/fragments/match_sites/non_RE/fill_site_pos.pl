use strict;
use warnings;

# scope:
#   this script is applied to non-RE-cleaved fragments, i.e., a non-agFree library,
#   to skip non-applicable actions in match_nodes_to_sites.pl
# action:
#   fill actual read endpoints instead of matching RE site positions
# input:
#   name-sorted SITE_SAM on STDIN
# output: 
#   name-sorted SITE_SAM on STDOUT
#   where this script updates values SITE_POS1_1 and SITE_POS1_2 as alignment node positions

# initialize reporting
our $script = "fill_site_pos";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);

# environment variables

# constants
use constant {
    S_QNAME             => 0, # SITE_SAM fields
    S_FLAG              => 1,
    S_RNAME             => 2,
    S_POS1              => 3, # 1-based
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    DE_TAG              => 6,
    CS_TAG              => 7,
    XF_TAG              => 8,
    XH_TAG              => 9,
    N_REF_BASES         => 10,
    N_READ_BASES        => 11,
    BLOCK_N             => 12,
    SITE_INDEX1_1       => 13,
    SITE_POS1_1         => 14,
    SITE_HAPS_1         => 15,
    SITE_DIST_1         => 16,
    SITE_INDEX1_2       => 17,
    SITE_POS1_2         => 18,
    SITE_HAPS_2         => 19,
    SITE_DIST_2         => 20,
    SEQ_SITE_INDEX1_2   => 21,
    SEQ_SITE_POS1_2     => 22,
    SEQ_SITE_HAPS_2     => 23,
    IS_END_TO_END       => 24,
    READ_HAS_JXN        => 25,
    TARGET_CLASS        => 26,
    S_SEQ               => 27,
    S_QUAL              => 28,
    #-------------
    SPLIT_TO_SITE_DIST_2 => 22,
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
};

# process alignments one at a time
while(my $aln = <STDIN>){

    # get closest RE sites
    # adjust alignment positions for strand/end offset but NOT clip yet as clip handling differs for outer and inner nodes
    # node1 as recorded here is always before node2 in read order
    #        *       *       *     proper sitePos1 values
    #   ----|5-----3/5-----3|----  as nodes are numbered for alignments
    #   ----|3-----5/3-----5|----
    my @aln = split("\t", $aln, SPLIT_TO_SITE_DIST_2);
    if($aln[S_FLAG] & _REVERSE){
        $aln[SITE_POS1_1] = getEnd(@aln[S_POS1, S_CIGAR]);
        $aln[SITE_POS1_2] = $aln[S_POS1]; # + 1
    } else {
        $aln[SITE_POS1_1] = $aln[S_POS1];
        $aln[SITE_POS1_2] = getEnd(@aln[S_POS1, S_CIGAR]); # + 1
    }
    print join("\t", @aln);
}

use strict;
use warnings;

# actions:
#   mask SEQ, QUAL and CIGAR from final junctions files
# input:
#   junctions on STDIN
# output:
#   junctions on STDOUT with same column format but SEQ, QUAL and CIGAR masked to *

# initialize reporting
our $script = "package";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);

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
    UJXN_READ_HAS_SV        => 12,
    UJXN_JSRC_MAPQ          => 13, 
    UJXN_JSRC_DE_TAG        => 14,
    UJXN_JSRC_SITE_DIST     => 15,
    UJXN_JSRC_QNAMES        => 16, 
    UJXN_JSRC_SEQS          => 17, # masked by this script
    UJXN_JSRC_QUALS         => 18, # masked by this script
    UJXN_JSRC_CIGARS        => 19, # masked by this script (included since ONT CIGAR can be quite long)
    UJXN_JSRC_ORIENTATIONS  => 20,
    UJXN_HAS_ALT_ALIGNMENT  => 21,
    UJXN_SV_SIZE            => 22,
    UJXN_IS_INTERGENOME     => 23,
    UJXN_TARGET_1           => 24,
    UJXN_TARGET_DIST_1      => 25,
    UJXN_TARGET_2           => 26,
    UJXN_TARGET_DIST_2      => 27,
    UJXN_GENES_1            => 28,
    UJXN_GENE_DIST_1        => 29,
    UJXN_GENES_2            => 30,
    UJXN_GENE_DIST_2        => 31,
    UJXN_IS_EXCLUDED_1      => 32,
    UJXN_IS_EXCLUDED_2      => 33,
    UJXN_BKPT_COVERAGE_1    => 34,
    UJXN_BKPT_COVERAGE_2    => 35,
    #-------------------
    SPLIT_TO_UJXN_JSRC_CIGARS => 21
};

# process junctions 
while (my $jxn = <STDIN>){
    my @jxn = split("\t", $jxn, SPLIT_TO_UJXN_JSRC_CIGARS);
    $jxn[UJXN_JSRC_SEQS] = '*';
    $jxn[UJXN_JSRC_QUALS] = '*';
    $jxn[UJXN_JSRC_CIGARS] = '*';
    print join("\t", @jxn);
}

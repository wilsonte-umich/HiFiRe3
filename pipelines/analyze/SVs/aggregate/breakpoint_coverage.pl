use strict;
use warnings;

# actions:
#   add information on the total read coverage crossing each breakpoint
# input:
#   output of junction_coverage.pl repeated as dictated by bedtools interect -wb on STDIN
# output:
#   unique junctions with all final metadata columns on STDOUT

# initialize reporting
our $script = "breakpoint_coverage";
our $error  = "$script error";
my @jxns;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
resetCountFile();

# environment variables
fillEnvVar(\our $DEDUPLICATE_READS, 'DEDUPLICATE_READS');

# constants
use constant {
    I_ALN_CHROM_INDEX1        => 0,
    I_ALN_START0              => 1,
    I_ALN_END1                => 2,
    I_ALN_N_OBSERVED          => 3,
    I_JXN_CHROM_INDEX1        => 4,
    I_JXN_BKPT_START0         => 5,
    I_JXN_BKPT_END1           => 6,
    I_JXN_JXN_I0              => 7,
    I_JXN_BKPT_I0             => 8,
    I_JXN_JXN                 => 9,
    #-------------
    I_SPLIT_TO_JXN_JXN        =>10, # leaves JXN concatenated at first split
    #-------------
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
    UJXN_JSRC_SEQS          => 17,
    UJXN_JSRC_QUALS         => 18,
    UJXN_JSRC_CIGARS        => 19,
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
    UJXN_BKPT_COVERAGE_1    => 34, # added by this script
    UJXN_BKPT_COVERAGE_2    => 35,
};

# process breakpoints one at a time
# expect multiple lines for the same breakpoint for every read that crossed it
# accumulate in memory while counting
while (my $intsct = <STDIN>){
    chomp $intsct;
    my @intsct = split("\t", $intsct, I_SPLIT_TO_JXN_JXN);
    unless($jxns[$intsct[I_JXN_JXN_I0]]){
        my @jxn = split("\t", $intsct[I_JXN_JXN]);
        $jxns[$intsct[I_JXN_JXN_I0]] = \@jxn;
    }
    my $j = $intsct[I_JXN_BKPT_I0] + UJXN_BKPT_COVERAGE_1;
    if ($DEDUPLICATE_READS){
        $jxns[$intsct[I_JXN_JXN_I0]][$j]++; # only partially deduplicated
    } else {
        $jxns[$intsct[I_JXN_JXN_I0]][$j] += $intsct[I_ALN_N_OBSERVED];
    }
}

# output unique junctions with all final metadata columns
for my $jxn(@jxns){
    $jxn or next;
    print join("\t", @$jxn), "\n";
}

# [wilsonte@gl-login2 tagFreeLR_pilot_WL002_082724]$ zcat site_sam/tagFreeLR_pilot_WL002_082724.GM12878_unspecified_EcoRV.chr1.site_sam.gz | grep 95a80e45-7504-4a3d-9aac-2c9fdde63e26:1028:0:0:0:0:0:0:1 | cut -f 1-5
# 95a80e45-7504-4a3d-9aac-2c9fdde63e26:1028:0:0:0:0:0:0:1 16      chr1    143701  60
# 95a80e45-7504-4a3d-9aac-2c9fdde63e26:1028:0:0:0:0:0:0:1 2048    chr13   104551419       60

# [wilsonte@gl-login2 tagFreeLR_pilot_WL002_082724]$ zcat tagFreeLR_pilot_WL002_082724.GM12878_unspecified_EcoRV.alignments.txt.bgz | grep -P "\t143701\t"
# 1       176825  143701  0       1       4       0       166     1       1       1       1

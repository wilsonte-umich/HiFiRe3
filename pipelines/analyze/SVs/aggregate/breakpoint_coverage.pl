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
fillEnvVar(\our $DATA_NAME, 'DATA_NAME');

# constants
use constant {
    I_ALN_CHROM_INDEX1        => 0, # from SV_ALIGNMENTS_FILE via bedtools intersect in aggregate_SVs.sh
    I_ALN_START0              => 1,
    I_ALN_END1                => 2,
    I_ALN_N_OBSERVED          => 3,
    I_JXN_CHROM_INDEX1        => 4, # junction metadata prepended to junction by junction_coverage.pl
    I_JXN_BKPT_START0         => 5,
    I_JXN_BKPT_END1           => 6,
    I_JXN_JXN_I0              => 7,
    I_JXN_BKPT_I0             => 8,
    I_JXN_JXN                 => 9, # unsplit full junction, columns below
    #-------------
    I_SPLIT_TO_JXN_JXN        =>10, # leaves JXN concatenated at first split
    #-------------
    OJXN_CHROM_INDEX1_1     => 0, # full junction reported by junction_coverage.pl
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
    UJXN_IS_EXPECTED          => 26,
    UJXN_HAS_ALT_ALIGNMENT    => 27,
    UJXN_SV_SIZE              => 28,
    UJXN_IS_INTERGENOME       => 29,
    UJXN_TARGET_1             => 30,
    UJXN_TARGET_DIST_1        => 31,
    UJXN_TARGET_2             => 32,
    UJXN_TARGET_DIST_2        => 33,
    UJXN_GENES_1              => 34,
    UJXN_GENE_DIST_1          => 35,
    UJXN_GENES_2              => 36,
    UJXN_GENE_DIST_2          => 37,
    UJXN_IS_EXCLUDED_1        => 38,
    UJXN_IS_EXCLUDED_2        => 39,
    UJXN_BKPT_COVERAGE_1      => 40, # added by this script
    UJXN_BKPT_COVERAGE_2      => 41,
    CMP_N_SAMPLES             => 42,
    CMP_SAMPLES               => 43,
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
    print join("\t", @$jxn, 1, ",$DATA_NAME,"), "\n";
}

use strict;
use warnings;

# actions:
#   collect additional sample SV statistics now that all filters and aggregation are done
# input:
#   unique junctions on STDIN
# output:
#   one line per final junction breakpoint on STDOUT, for use in breakpoint coverage assessment
#   statistics reported to log
#   ANALYSIS_CHROMS_FILE, to support app conversion between chrom and chromIndex

# initialize reporting
our $script = "junction_coverage";
our $error  = "$script error";
my ($nInputJxns, $nExpectedJxns) = (0) x 10;
my %jxnCoverage;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms targets genes exclude);
resetCountFile();

# environment variables
fillEnvVar(\our $ANALYSIS_CHROMS_FILE, 'ANALYSIS_CHROMS_FILE');
fillEnvVar(\our $GENOME_FASTA,         'GENOME_FASTA');
fillEnvVar(\our $IS_COMPOSITE_GENOME,  'IS_COMPOSITE_GENOME');
fillEnvVar(\our $TARGETS_BED,          'TARGETS_BED');
fillEnvVar(\our $REGION_PADDING,       'REGION_PADDING');
fillEnvVar(\our $GENES_BED,            'GENES_BED');
fillEnvVar(\our $GENOME_EXCLUSIONS_BED,'GENOME_EXCLUSIONS_BED');
our $TARGET_SCALAR = 50;
our $GENE_SCALAR = 50;

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
    UJXN_IS_EXPECTED          => 26,
    UJXN_HAS_ALT_ALIGNMENT    => 27,
    UJXN_SV_SIZE              => 28, # added by this script
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
    UJXN_BKPT_COVERAGE_1      => 40, # added downstream, not present yet
    UJXN_BKPT_COVERAGE_2      => 41,
    CMP_N_SAMPLES             => 42,
    CMP_SAMPLES               => 43,
    #-------------
    EXPECTED   => 1, # junction sequencing orientations
    UNEXPECTED => 0,
    #-------------
    TYPE_PROPER         => 0, # junction/edge types
    TYPE_DELETION       => 1,
    TYPE_DUPLICATION    => 2,
    TYPE_INVERSION      => 3,
    TYPE_TRANSLOCATION  => 4,
    #-------------
    TRUE  => 1,
    FALSE => 0,
};

# initialize the genome
use vars qw(%revChromIndex);
setCanonicalChroms();

# initialize the target regions and genes, if any
use vars qw($nRegions $nGenes);
loadTargetRegions();
loadGeneRegions();
initializeExclude($GENOME_EXCLUSIONS_BED);

# write the chroms file
writeChromsFile($ANALYSIS_CHROMS_FILE, $GENOME_FASTA);

# process junctions one key at a time
while (my $jxn = <STDIN>){
    chomp $jxn;
    my @jxn = split("\t", $jxn);
    $nInputJxns++;
    $jxn[UJXN_IS_EXPECTED] and $nExpectedJxns++;
    $jxnCoverage{$jxn[UJXN_N_OBSERVED]}[$jxn[UJXN_IS_EXPECTED]]++;

    $jxn[UJXN_SV_SIZE] = 
        $jxn[OJXN_JXN_TYPE] == TYPE_TRANSLOCATION ? 
        0 :
        abs($jxn[OJXN_REF_POS1_2] - $jxn[OJXN_REF_POS1_1]);

    $jxn[UJXN_IS_INTERGENOME] = (
        $IS_COMPOSITE_GENOME and 
        $jxn[OJXN_JXN_TYPE] == TYPE_TRANSLOCATION and 
        isInterGenome(@jxn[OJXN_CHROM_INDEX1_1, OJXN_CHROM_INDEX1_2])
    ) ? TRUE : FALSE;

    ($jxn[UJXN_TARGET_1], $jxn[UJXN_TARGET_DIST_1]) = 
        getPosTarget($revChromIndex{$jxn[OJXN_CHROM_INDEX1_1]}, $jxn[OJXN_REF_POS1_1]);
    ($jxn[UJXN_TARGET_2], $jxn[UJXN_TARGET_DIST_2]) = 
        getPosTarget($revChromIndex{$jxn[OJXN_CHROM_INDEX1_2]}, $jxn[OJXN_REF_POS1_2]);
    ($jxn[UJXN_GENES_1], $jxn[UJXN_GENE_DIST_1]) = 
        getPosGene($revChromIndex{$jxn[OJXN_CHROM_INDEX1_1]}, $jxn[OJXN_REF_POS1_1]);
    ($jxn[UJXN_GENES_2], $jxn[UJXN_GENE_DIST_2]) = 
        getPosGene($revChromIndex{$jxn[OJXN_CHROM_INDEX1_2]}, $jxn[OJXN_REF_POS1_2]);

    $jxn[UJXN_IS_EXCLUDED_1] = 
        isExcludedPosition($jxn[OJXN_CHROM_INDEX1_1], $jxn[OJXN_REF_POS1_1]);
    $jxn[UJXN_IS_EXCLUDED_2] = 
        isExcludedPosition($jxn[OJXN_CHROM_INDEX1_2], $jxn[OJXN_REF_POS1_2]);

    # repeat junction for use for breakpoint coverage determination, one line per breakpoint
    my $jxn = join("\t", @jxn)."\n";
    print join("\t", $jxn[OJXN_CHROM_INDEX1_1], $jxn[OJXN_REF_POS1_1] - 1, $jxn[OJXN_REF_POS1_1], $nInputJxns - 1, 0, $jxn);
    print join("\t", $jxn[OJXN_CHROM_INDEX1_2], $jxn[OJXN_REF_POS1_2] - 1, $jxn[OJXN_REF_POS1_2], $nInputJxns - 1, 1, $jxn);
}

# print summary information
printCount(commify($nInputJxns),    'nInputJxns',    'input junctions');
# printCount(commify($nExpectedJxns), 'nExpectedJxns', 'junctions found in at least one read where all SVs resolved in a haplotype alignment');
print STDERR join("\t", qw(jxnCoverage nInExpectedPaths nInUnexpectedPaths)), "\n";
foreach my $jxnCoverage(sort { $a <=> $b} keys %jxnCoverage){
    print STDERR join("\t", 
        $jxnCoverage, 
        $jxnCoverage{$jxnCoverage}[  EXPECTED] || 0, 
        $jxnCoverage{$jxnCoverage}[UNEXPECTED] || 0
    ), "\n"
}

sub isInterGenome {
    my ($chromIndex1_1, $chromIndex1_2) = @_;
    my ($genome1) = ($revChromIndex{$chromIndex1_1} =~ m/.+_(.+)/); # e.g., chr1_(hs1)
    my ($genome2) = ($revChromIndex{$chromIndex1_2} =~ m/.+_(.+)/);
    return $genome1 ne $genome2;
}

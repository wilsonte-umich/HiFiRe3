use strict;
use warnings;

# actions:
#   purge PCR duplicates when junctions appeared on sheared (not RE-cleaved) molecules with the same outer endpoints
#   applies to tagFree libraries with randomly sheared fragments, not ligFree libraries with fixed fragment endpoints
# input:
#   unique junctions with merged source metadata on STDIN
# output:
#   unique junctions with metadata adjusted for PCR duplicates on STDOUT

# initialize reporting
our $script = "purge_PCR duplicates";
our $error  = "$script error";
my ($nInputJxns, $nInputEvents, $nOutputEvents) = (0) x 10;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
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
    UJXN_READ_HAS_SV        => 12,
    UJXN_JSRC_MAPQ          => 13, 
    UJXN_JSRC_DE_TAG        => 14,
    UJXN_JSRC_SITE_DIST     => 15,
    UJXN_JSRC_QNAMES        => 16, 
    UJXN_JSRC_SEQS          => 17,
    UJXN_JSRC_QUALS         => 18,
    UJXN_JSRC_CIGARS        => 19,
    UJXN_JSRC_ORIENTATIONS  => 20,
    UJXN_JSRC_CHANNELS      => 21, # used and then dropped by this script
    UJXN_JSRC_OUTER_ENDPOINTS => 22,
};

# process junctions one key at a time
while (my $jxn = <STDIN>){
    chomp $jxn;
    my @jxn = split("\t", $jxn);
    $nInputJxns++;
    $nInputEvents += $jxn[UJXN_N_OBSERVED];

    # TODO: if needed, could implement fuzzing matching of outer endpoint positions
    my @outerOEs = split(",", $jxn[UJXN_JSRC_OUTER_ENDPOINTS]);
    my %outerOEs = map { $_ => 1 } @outerOEs;

    $jxn[UJXN_N_OBSERVED] = scalar(keys %outerOEs);
    $nOutputEvents += $jxn[UJXN_N_OBSERVED];
    print join("\t", @jxn[OJXN_CHROM_INDEX1_1..UJXN_JSRC_ORIENTATIONS]), "\n";
}

# print summary information
printCount(commify($nInputJxns),    'nInputJxns',    'input junctions');
printCount(commify($nInputEvents),  'nInputEvents',  'total sequencing events supporting nInputJxns');
printCount(commify($nOutputEvents), 'nOutputEvents', 'sequencing events supporting nInputJxns after discarding PCR duplicates');

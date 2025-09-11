use strict;
use warnings;

# actions:
#   purge ONT duplexes when junctions appeared on different strands in the same channel
#   junctions sequenced on the same strand in the same channel, or in different channels
#       must be independent read sequencing events in a non-amplified library
# input:
#   unique junctions with merged source metadata on STDIN
# output:
#   unique junctions with metadata adjusted for ONT duplex duplicates on STDOUT

# initialize reporting
our $script = "purge_ONT_duplexes";
our $error  = "$script error";
my ($nInputJxns, $nDuplexes, $nInputEvents, $nOutputEvents) = (0) x 10;
my %duplexJxns;

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
    #-----------------
    CANONICAL    => 0, # junction sequencing orientations
    NONCANONICAL => 1
};

# process junctions one key at a time
while (my $jxn = <STDIN>){
    chomp $jxn;
    my @jxn = split("\t", $jxn);
    $nInputJxns++;
    $nInputEvents += $jxn[UJXN_N_OBSERVED];
    my @orientations = split(",", $jxn[UJXN_JSRC_ORIENTATIONS]);
    my @channels     = split(",", $jxn[UJXN_JSRC_CHANNELS]);
    my %channelOrientations;
    foreach my $i(0..$#channels){
        $channelOrientations{$channels[$i]}{$orientations[$i]}++;
    }
    $jxn[UJXN_N_OBSERVED] = 0; 
    my %jxnOrientations;
    foreach my $channel(keys %channelOrientations){
        my @channelOrientations = keys %{$channelOrientations{$channel}};
        if(@channelOrientations > 1){
            my @oris = sort { $channelOrientations{$channel}{$b} <=> $channelOrientations{$channel}{$a} } CANONICAL, NONCANONICAL;
            $jxn[UJXN_N_OBSERVED] += $channelOrientations{$channel}{$oris[0]};
            $jxnOrientations{$oris[0]}++;
            $duplexJxns{$jxn}++;
            $nDuplexes++;
        } else {
            $jxn[UJXN_N_OBSERVED] += $channelOrientations{$channel}{$channelOrientations[0]};
            $jxnOrientations{$channelOrientations[0]}++;
        }
    }
    $nOutputEvents += $jxn[UJXN_N_OBSERVED];
    print join("\t", @jxn[OJXN_CHROM_INDEX1_1..UJXN_JSRC_ORIENTATIONS]), "\n";
}

# print summary information
printCount(commify($nInputJxns),    'nInputJxns',    'input junctions');
printCount(commify(scalar(keys %duplexJxns)), 'nDuplexJxns',   'input junctions with at least one ONT duplex event pair');
printCount(commify($nDuplexes),     'nDuplexes',     'duplex events with two-stranded junction detection in one channel');
printCount(commify($nInputEvents),  'nInputEvents',  'total sequencing events supporting nInputJxns');
printCount(commify($nOutputEvents), 'nOutputEvents', 'sequencing events supporting nInputJxns after discarding duplex duplicates');

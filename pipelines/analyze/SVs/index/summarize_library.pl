use strict;
use warnings;

# actions:
#   calculate read and base counts over a (deduplicated) library
# input:
#   chromIndex1_1,chromIndex1_2,(reordered)alnI0,sumRefBases,sumReadBases,nAlns on STDIN
# output:
#   various counts in a log report

# initialize reporting
our $script = "summarize_library";
our $error  = "$script error";
my ($nEvents, $nAlns, $nRefBases, $nReadBases) = (0) x 100;
my (%nEvents, %nRefBases, %nReadBases);

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
resetCountFile();

# environment variables
fillEnvVar(\our $GENOME,              'GENOME');
fillEnvVar(\our $IS_COMPOSITE_GENOME, 'IS_COMPOSITE_GENOME');

# constants
use constant {
    CHROM_INDEX1_1  => 0,
    CHROM_INDEX1_2  => 1,
    ALN_I0          => 2,
    SUM_REF_BASES   => 3,
    SUM_READ_BASES  => 4,
    N_GROUP_ALNS    => 5,
};

# initialize the genome
use vars qw(%revChromIndex);
setCanonicalChroms();

# tally reads and bases overall and per composite genome
while (my $grp = <STDIN>){
    chomp $grp;
    my @grp = split("\t", $grp);
    $grp[ALN_I0] == 0 and $nEvents += $grp[N_GROUP_ALNS];
    $nAlns      += $grp[N_GROUP_ALNS];
    $nRefBases  += $grp[SUM_REF_BASES];
    $nReadBases += $grp[SUM_READ_BASES];

    # count reads and bases per genome to determine observed mixing ratio
    if($IS_COMPOSITE_GENOME and $grp[CHROM_INDEX1_1] == $grp[CHROM_INDEX1_2]){
        my ($genome) = ($revChromIndex{$grp[CHROM_INDEX1_1]} =~ m/.+_(.+)/); # e.g., chr1_(hs1)
        $grp[ALN_I0] == 0 and $nEvents{$genome} += $grp[N_GROUP_ALNS];
        $nRefBases{$genome}  += $grp[SUM_REF_BASES];
        $nReadBases{$genome} += $grp[SUM_READ_BASES];
    }
}

# report counts
printCount(commify($nEvents),       'nEvents',      'final (deduplicated) read count');
printCount(commify($nAlns),         'nAlns',        'final (deduplicated) alignment count');
printCount(commify($nReadBases),    'nReadBases',   'final (deduplicated) read bases in alignments');
printCount(commify($nRefBases),     'nRefBases',    'final (deduplicated) reference bases in alignments');
if($IS_COMPOSITE_GENOME){
    for my $genome (keys %nEvents){
        printCount(commify($nEvents{$genome}),       "nEvents_$genome",      "final (deduplicated) read count, $genome (w/o translocations)");
        printCount(commify($nReadBases{$genome}),    "nReadBases_$genome",   "final (deduplicated) read bases in alignments, $genome (w/o translocations)");
        printCount(commify($nRefBases{$genome}),     "nRefBases_$genome",    "final (deduplicated) reference bases in alignments, $genome (w/o translocations)");
    }
}

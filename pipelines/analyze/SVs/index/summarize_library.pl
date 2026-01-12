use strict;
use warnings;

# actions:
#   calculate read and base counts over a (deduplicated) library
#   only on-target reads are counted since off-target reads were filtered upstream by analyze/fragments/apply_filters.pl
# input:
#   orderedNode1,orderedNode2,channel,(reordered)alnI0,nRefBases,nReadBases on STDIN
# output:
#   various counts in a log report

# initialize reporting
our $script = "summarize_library";
our $error  = "$script error";
my ($nReadsOn,  $nAlnsOn,  $nRefBasesOn,  $nReadBasesOn) = (0) x 100;
my (%nReadsOn,  %nRefBasesOn,  %nReadBasesOn);

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
    ORDERED_NODE_1 => 0,
    ORDERED_NODE_2 => 1,
    CHANNEL        => 2,
    ALN_I0         => 3,
    N_REF_BASES    => 4,
    N_READ_BASES   => 5,
};

# initialize the genome
use vars qw(%revChromIndex);
setCanonicalChroms();

# tally reads and bases overall and per composite genome
while (my $grp = <STDIN>){
    chomp $grp;
    my @grp = split("\t", $grp);
    $grp[ALN_I0] == 0 and $nReadsOn++;
    $nAlnsOn++;
    $nRefBasesOn  += $grp[N_REF_BASES];
    $nReadBasesOn += $grp[N_READ_BASES];

    # count reads and bases per genome to determine observed mixing ratio
    if($IS_COMPOSITE_GENOME){
        my $chrom1 = chromFromNode($grp[ORDERED_NODE_1]);
        my $chrom2 = chromFromNode($grp[ORDERED_NODE_2]);
        if($chrom1 eq $chrom2){
            my ($genome) = ($chrom1 =~ m/.+_(.+)/); # e.g., chr1_(hs1)
            $grp[ALN_I0] == 0 and $nReadsOn{$genome}++;
            $nRefBasesOn{$genome}  += $grp[N_REF_BASES];
            $nReadBasesOn{$genome} += $grp[N_READ_BASES];
        } 
    }
}
sub chromFromNode {
    my ($node) = @_;
    $node =~ m/(\d{1,2})\d{9}$/;
    $revChromIndex{$1}
}

# report counts
printCount(commify($nReadsOn),        'nReadsOn',       '(deduped) read count, 5\' on target (or whole genome)');
printCount(commify($nAlnsOn),         'nAlnsOn',        '(deduped) alignment count, 5\' on target (or whole genome)');
printCount(commify($nReadBasesOn),    'nReadBasesOn',   '(deduped) read bases in alignments, 5\' on target (or whole genome)');
printCount(commify($nRefBasesOn),     'nRefBasesOn',    '(deduped) reference bases in alignments, 5\' on target (or whole genome)');
if($IS_COMPOSITE_GENOME){
    for my $genome (keys %nReadsOn){
        printCount(commify($nReadsOn{$genome}),        "nReadsOn_$genome",       "(deduped) read count, $genome (w/o transloc), 5\' on target");
        printCount(commify($nReadBasesOn{$genome}),    "nReadBasesOn_$genome",   "(deduped) read bases in alignments, $genome (w/o transloc), 5\' on target");
        printCount(commify($nRefBasesOn{$genome}),     "nRefBasesOn_$genome",    "(deduped) reference bases in alignments, $genome (w/o transloc), 5\' on target");
    }
}

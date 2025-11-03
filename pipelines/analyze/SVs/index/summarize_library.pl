use strict;
use warnings;

# actions:
#   calculate read and base counts over a (deduplicated) library
# input:
#   chromIndex1_1,chromIndex1_2,end5OnTarget,(reordered)alnI0,sumRefBases,sumReadBases,nAlns on STDIN
# output:
#   various counts in a log report

# initialize reporting
our $script = "summarize_library";
our $error  = "$script error";
my (
    $nEventsOn,  $nAlnsOn,  $nRefBasesOn,  $nReadBasesOn,
    # $nEventsOff, $nAlnsOff, $nRefBasesOff, $nReadBasesOff
) = (0) x 100;
my (
    %nEventsOn,  %nRefBasesOn,  %nReadBasesOn,
    # %nEventsOff, %nRefBasesOff, %nReadBasesOff
);

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
    # END5_ON_TARGET  => 2, # all indexed reads are now on-target for targeted libraries
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
    # if($grp[END5_ON_TARGET]){
        $grp[ALN_I0] == 0 and $nEventsOn += $grp[N_GROUP_ALNS];
        $nAlnsOn      += $grp[N_GROUP_ALNS];
        $nRefBasesOn  += $grp[SUM_REF_BASES];
        $nReadBasesOn += $grp[SUM_READ_BASES];
    # } else {
        # $grp[ALN_I0] == 0 and $nEventsOff += $grp[N_GROUP_ALNS];
        # $nAlnsOff      += $grp[N_GROUP_ALNS];
        # $nRefBasesOff  += $grp[SUM_REF_BASES];
        # $nReadBasesOff += $grp[SUM_READ_BASES];
    # }

    # count reads and bases per genome to determine observed mixing ratio
    if($IS_COMPOSITE_GENOME and $grp[CHROM_INDEX1_1] == $grp[CHROM_INDEX1_2]){
        my ($genome) = ($revChromIndex{$grp[CHROM_INDEX1_1]} =~ m/.+_(.+)/); # e.g., chr1_(hs1)
        # if($grp[END5_ON_TARGET]){
            $grp[ALN_I0] == 0 and $nEventsOn{$genome} += $grp[N_GROUP_ALNS];
            $nRefBasesOn{$genome}  += $grp[SUM_REF_BASES];
            $nReadBasesOn{$genome} += $grp[SUM_READ_BASES];
        # } else {
        #     $grp[ALN_I0] == 0 and $nEventsOff{$genome} += $grp[N_GROUP_ALNS];
        #     $nRefBasesOff{$genome}  += $grp[SUM_REF_BASES];
        #     $nReadBasesOff{$genome} += $grp[SUM_READ_BASES];
        # }
    }
}

# report counts
printCount(commify($nEventsOn),       'nEventsOn',      '(deduped) read count, 5\' on target (or whole genome)');
printCount(commify($nAlnsOn),         'nAlnsOn',        '(deduped) alignment count, 5\' on target (or whole genome)');
printCount(commify($nReadBasesOn),    'nReadBasesOn',   '(deduped) read bases in alignments, 5\' on target (or whole genome)');
printCount(commify($nRefBasesOn),     'nRefBasesOn',    '(deduped) reference bases in alignments, 5\' on target (or whole genome)');
# printCount(commify($nEventsOff),      'nEventsOff',     '(deduped) read count, 5\' off target');
# printCount(commify($nAlnsOff),        'nAlnsOff',       '(deduped) alignment count, 5\' off target');
# printCount(commify($nReadBasesOff),   'nReadBasesOff',  '(deduped) read bases in alignments, 5\' off target');
# printCount(commify($nRefBasesOff),    'nRefBasesOff',   '(deduped) reference bases in alignments, 5\' off target');
if($IS_COMPOSITE_GENOME){
    for my $genome (keys %nEventsOn){
        printCount(commify($nEventsOn{$genome}),       "nEventsOn_$genome",      "(deduped) read count, $genome (w/o transloc), 5\' on target");
        printCount(commify($nReadBasesOn{$genome}),    "nReadBasesOn_$genome",   "(deduped) read bases in alignments, $genome (w/o transloc), 5\' on target");
        printCount(commify($nRefBasesOn{$genome}),     "nRefBasesOn_$genome",    "(deduped) reference bases in alignments, $genome (w/o transloc), 5\' on target");
    }
    # for my $genome (keys %nEventsOff){
    #     printCount(commify($nEventsOff{$genome}),       "nEventsOff_$genome",      "(deduped) read count, $genome (w/o transloc), 5\' off target");
    #     printCount(commify($nReadBasesOff{$genome}),    "nReadBasesOff_$genome",   "(deduped) read bases in alignments, $genome (w/o transloc), 5\' off target");
    #     printCount(commify($nRefBasesOff{$genome}),     "nRefBasesOff_$genome",    "(deduped) reference bases in alignments, $genome (w/o transloc), 5\' off target");
    # }
}

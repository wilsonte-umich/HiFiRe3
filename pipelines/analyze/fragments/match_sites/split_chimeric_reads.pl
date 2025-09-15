use strict;
use warnings;

# action:
#   split sequences/molecules based on node distance to closest RE site
#       junctions are rejected if either node is <= REJECT_JUNCTION_DISTANCE bp from an RE site
#       the flanking alignments are consistent with a foldback inversion 
#       adapter is present in junction inserted bases of sufficient length
#   if ONT, additionally split sequences/molecules to contigs if
#       junction inserted bases have a low quality span, e.g., as results from an ONT follow-on event
#   when a junction is rejected as chimeric, discard all alignments more 3'/distal to that junction (and thus the junction itself)
# input:
#   name-sorted SITE_SAM on STDIN
# output: 
#   name-sorted SITE_SAM on STDOUT
#   where this script may remove 3'/distal read alignments but does not modify those that remain

# initialize reporting
our $script = "split_chimeric_reads";
our $error  = "$script error";
my ($nKeptJxns, $nRejectedJxns, $nFoldbackRejections, $nInsertionCuts, $nAdapterRejections) = (0) x 10;
my (@siteDistances, @adapterScores);

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general smith_waterman);
resetCountFile();

# environment variables
fillEnvVar(\our $SEQUENCING_PLATFORM,       'SEQUENCING_PLATFORM');
fillEnvVar(\our $REJECT_JUNCTION_DISTANCE,  'REJECT_JUNCTION_DISTANCE');
fillEnvVar(\our $INSERTION_WINDOW_SIZE,     'INSERTION_WINDOW_SIZE');
fillEnvVar(\our $MIN_INSERTION_WINDOW_QUAL, 'MIN_INSERTION_WINDOW_QUAL');
fillEnvVar(\our $INSERTION_ADAPTER_SEQUENCE,'INSERTION_ADAPTER_SEQUENCE');
fillEnvVar(\our $MIN_ADAPTER_LENGTH,        'MIN_ADAPTER_LENGTH'); # only look for adapters in junction inserts in this size range
fillEnvVar(\our $MAX_ADAPTER_LENGTH,        'MAX_ADAPTER_LENGTH');
fillEnvVar(\our $MIN_ADAPTER_SCORE,         'MIN_ADAPTER_SCORE');  # empirically determined SW score that signals presence of an adapter in a junction insertion
my $minSumInsQual = ($MIN_INSERTION_WINDOW_QUAL + 33) * $INSERTION_WINDOW_SIZE;

# initialize the genome files
use vars qw(%chromIndex @canonicalChroms);
setCanonicalChroms();

# set platform-specific parameters
my $isONT = $SEQUENCING_PLATFORM eq "ONT";
my $adapterCore = $INSERTION_ADAPTER_SEQUENCE; # duplex portion of the ONT kit adapter; for ligation kit, last T matches the one-base A-tail; fused to 5' genomic ends 
my $adapterCoreRc = $adapterCore;    # fused to 3' genomic ends
rc(\$adapterCoreRc);

# constants
use constant {
    S_QNAME             => 0, # SITE_SAM fields
    S_FLAG              => 1,
    S_RNAME             => 2,
    S_POS1              => 3, # 1-based
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    CH_TAG              => 6,
    TL_TAG              => 7,
    DE_TAG              => 8,
    HV_TAG              => 9,
    N_REF_BASES         => 10,
    N_READ_BASES        => 11,
    BLOCK_N             => 12,
    SITE_INDEX1_1       => 13,
    SITE_POS1_1         => 14,
    SITE_DIST_1         => 15,
    SITE_INDEX1_2       => 16,
    SITE_POS1_2         => 17,
    SITE_DIST_2         => 18,
    SEQ_SITE_INDEX1_2   => 19,
    SEQ_SITE_POS1_2     => 20,
    IS_END_TO_END       => 21,
    READ_HAS_JXN        => 22,
    TARGET_CLASS        => 23,
    S_SEQ               => 24,
    S_QUAL              => 25,
    CS_TAG              => 26,
    #-------------
    SPLIT_TO_S_QUAL     => 27,
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
    #--------------
    MAX_SITE_DISTANCE => 100, # just used for tabulating junction site distances, not a filtering threshold
    #--------------
    FALSE => 0,
    TRUE  => 1
};

# process alignments one read at a time
my ($prevRead, @alns);
my ($clip1, $clip2, $insSize);
while(my $aln = <STDIN>){
    my @aln = split("\t", $aln, SPLIT_TO_S_QUAL);
    if($prevRead and $prevRead ne $aln[S_QNAME]){
        processRead();
        @alns = ();
    }
    push @alns, \@aln;
    $prevRead = $aln[S_QNAME];
}
processRead();

# process sequences by QNAME
sub processRead {

    # first, 5'-most read alignment on a canonical chrom is always valid (as checked upstream) and passed as is
    $chromIndex{$alns[0][S_RNAME]} or return;
    print join("\t", @{$alns[0]});

    # only SV-containing reads have additional junction:alignment(s) to process
    @alns > 1 or return;
    foreach my $i(1..$#alns){

        # reject junctions that are not on a canonical chrom
        if(!$chromIndex{$alns[$i][S_RNAME]}){
            $nRejectedJxns++;
            return; # stop checking further alignment in case read has more junctions to come
        }

        # calculate junction node distances
        # unlike outermost end matching, do not use clip-adjustment here (clips are large and reflect the SV)
        my $dist1 = abs( $alns[$i - 1][SITE_DIST_2] ); # junction node 1, at query 3' end of left  flanking alignment
        my $dist2 = abs( $alns[$i    ][SITE_DIST_1] ); # junction node 2, at query 5' end of right flanking alignment
        $siteDistances[min($dist1, MAX_SITE_DISTANCE)]++;
        $siteDistances[min($dist2, MAX_SITE_DISTANCE)]++;

        # reject junctions that fail REJECT_JUNCTION_DISTANCE
        # reject ONT junctions that contain adapters, if not already caught by REJECT_JUNCTION_DISTANCE
        if(
            $dist1 <= $REJECT_JUNCTION_DISTANCE or 
            $dist2 <= $REJECT_JUNCTION_DISTANCE or
            isFoldbackInversion($alns[$i - 1], $alns[$i]) or
            isOntFollowOn($alns[$i - 1], $alns[$i]) or
            junctionHasAdapters($alns[$i - 1], $alns[$i])
        ){
            $nRejectedJxns++;
            return; # stop checking further alignment in case read has more junctions to come

        # keep junctions that pass REJECT_JUNCTION_DISTANCE and ONT adapter check
        } else {
            $nKeptJxns++;
            print join("\t", @{$alns[$i]});
        }
    }
}

# print summary information
printCount(commify($nKeptJxns),          'nKeptJxns',            'junctions not found to be chimeric');
printCount(commify($nRejectedJxns),      'nRejectedJxns',        'junctions found to be chimeric');
printCount(commify($nFoldbackRejections),'nFoldbackRejections',  'portion of nRejectedJxns rejected due to foldback inversion');
printCount(commify($nInsertionCuts),     'nInsertionCuts',       'portion of nRejectedJxns rejected due to stretch of low-quality inserted bases');
printCount(commify($nAdapterRejections), 'nAdapterRejections',   'portion of nRejectedJxns rejected due to ONT adapter in junction');

# print site distances to log
print STDERR "\njunction node-to-site distances\n";
foreach my $dist(0..MAX_SITE_DISTANCE){
    print STDERR join("\t", $dist, $siteDistances[$dist] || 0), "\n";
}
# print adapter SW scores to log
print STDERR "\nadapter SW scores\n";
foreach my $score(0..100){
    print STDERR join("\t", $score, $adapterScores[$score] || 0), "\n";
}

# check whether a sequence is consistent with an ONT duplex read as a single foldback inversion
# it is most sensitive and acceptable to reject any inversion junction with reverse-complement overlap between its flanking alignments
# UPDATE: do this for all libraries now
#   ----->
#         | inversion junction
#   <-----
sub isFoldbackInversion {
    my ($aln1, $aln2) = @_;
    $$aln1[S_RNAME] eq $$aln2[S_RNAME] or return FALSE; # translocation
    ($$aln1[S_FLAG] & _REVERSE) != ($$aln2[S_FLAG] & _REVERSE) or return FALSE; # deletion or duplication
    (
        $$aln1[S_POS1] <= getEnd(@{$aln2}[S_POS1, S_CIGAR]) and
        $$aln2[S_POS1] <= getEnd(@{$aln1}[S_POS1, S_CIGAR])
    ) or return FALSE;
    $nFoldbackRejections++;
    TRUE;
}

# determine if an ONT junction has a very low quality insertion span, identifying it as a two-insert follow-on event
sub isOntFollowOn {
    my ($aln1, $aln2) = @_; # alignments to the 5'/left and 3'/right of a junction, respectively

    # locate the insertion
    # stop if insertion of insufficient size to warrant adapter/follow-on detection
    # do this in advance of junctionHasAdapters even if not an ONT library
    $clip1 = 
        ($$aln1[S_FLAG] & _REVERSE) ? 
        getLeftClip($$aln1[S_CIGAR]) :
        getRightClip($$aln1[S_CIGAR]);
    $clip2 = 
        ($$aln2[S_FLAG] & _REVERSE) ? 
        getRightClip($$aln2[S_CIGAR]) :
        getLeftClip($$aln2[S_CIGAR]);
    $insSize = ($clip1 + $clip2) - length($$aln1[S_SEQ]);
    $isONT or return FALSE;
    $insSize < $MIN_ADAPTER_LENGTH and return FALSE;
    $insSize > $INSERTION_WINDOW_SIZE or return FALSE;

    # examine the inserted bases for low-quality base stretches as occurs during follow-on missing signals
    my $insertedQuals = 
        ($$aln1[S_FLAG] & _REVERSE) ? 
        substr(substr($$aln1[S_QUAL], 0, $clip1), -$insSize) :
        substr(substr($$aln1[S_QUAL], -$clip1), 0, $insSize);
    my @insertedQuals = split("", $insertedQuals);
    my $sumInsQual = 0;
    map{ $sumInsQual += ord($_) } @insertedQuals[0..($INSERTION_WINDOW_SIZE - 1)];
    if($sumInsQual < $minSumInsQual){
        $nInsertionCuts++;
        return TRUE;
    }
    # examine windows of bases to find no-base signal stretches
    # observed behavior is read1 ... (3' adapter) ... no-base signals/bases ... 5' adapter ... read2
    foreach my $i(1..($insSize - $INSERTION_WINDOW_SIZE)){
        $sumInsQual -= ord($insertedQuals[$i - 1]);
        $sumInsQual += ord($insertedQuals[$i + $INSERTION_WINDOW_SIZE - 1]);
        if($sumInsQual < $minSumInsQual){
            $nInsertionCuts++;
            return TRUE;
        }
    }
    return FALSE;
}

# determine if an ONT junction has adapters, identifying it as a two-insert event
sub junctionHasAdapters {
    my ($aln1, $aln2) = @_; # alignments to the 5'/left and 3'/right of a junction, respectively
    $insSize < $MIN_ADAPTER_LENGTH and return FALSE;
    $insSize > $MAX_ADAPTER_LENGTH and return FALSE;

    # next examine the inserted bases for ONT adapters using Smith-Waterman on both strands
    my $insertedBases = 
        ($$aln1[S_FLAG] & _REVERSE) ? 
        substr(substr($$aln1[S_SEQ], 0, $clip1), -$insSize) :
        substr(substr($$aln1[S_SEQ], -$clip1), 0, $insSize);
    my ($qryOnRef1, $score1) = smith_waterman($insertedBases, $adapterCore);
    $adapterScores[max(0, min($score1, 100))]++;
    if($score1 >= $MIN_ADAPTER_SCORE){
        $nAdapterRejections++;
        return TRUE;
    }
    my ($qryOnRef2, $score2) = smith_waterman($insertedBases, $adapterCoreRc);
    $adapterScores[max(0, min($score2, 100))]++;
    if($score2 >= $MIN_ADAPTER_SCORE){
        $nAdapterRejections++;
        return TRUE;
    }
    return FALSE;
}

# ONT adapter SW scores
# 0       270
# 1       441
# 2       138
# 3       39
# 4       27
# 5       13
# 6       7
# 7       15
# 8       8
# 9       9
# 10      6
# 11      6
#-----------
# 12      5
# 13      7
# 14      8
# 15      6
# 16      7
# 17      11
# 18      19
# 19      30
# 20      28
# 21      10
# 22      34
# 23      67

# 0       147
# 1       297
# 2       75
# 3       18
# 4       7
# 5       6
# 6       7
# 7       1
# 8       8
# 9       1
# 10      2
# 11      1
#-----------
# 12      1
# 13      3
# 14      3
# 15      0
# 16      4
# 17      0
# 18      8
# 19      12
# 20      10
# 21      2
# 22      11
# 23      41

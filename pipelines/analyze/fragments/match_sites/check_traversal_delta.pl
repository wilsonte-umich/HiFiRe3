use strict;
use warnings;

# action:
#   calculate the reference vs. read traversal delta between all sets of alignments (across all sets of junctions)
#   use calculated deltas and MIN_TRAVERSAL_DELTA to set BLOCK_N of all alignments
#       any set of alignments where traversalDelta < MIN_TRAVERSAL_DELTA are assigned the same BLOCK_N
#   junctions between alignments with the same BLOCK_N may be considered artifactual during SV analysis
#   block numbering is done before alignment quality truncation when all read alignments are still present
#       unlike other SV metrics, block numbering can depend on alignments distant from the two flanking a junction if nAln > 2
#       we want to ensure the even lower quality alignments have an opporunity to question a read's traversal delta
#   traversal deltas will fail and suppress false junctions when:
#       a large, low quality base string is called as a deletion SV with an equally large insertion
#       a similar base string is falsely aligned to another genomic location when it should have aligned inline
#           =========*==========*========= reference locus 1
#           ||||||||||          ||||||||||
#           ---------*~~~~~~~~~~*--------- query read, where lower quality ~~ bases are either an unaligned insertion or improperly aligned 
#                     ||||||||||
#           ============================== reference locus 2
#       notice that the traversal on reference and query are the same from * to *
#       contrast the above with a true deletion with some smaller number of bases inserted at the junction
#           =========*==========*========= reference locus 1
#           ||||||||||          ||||||||||
#           ---------*    ~~    *--------- query read, where lower quality ~~ bases are either an unmapped insertion or improperly aligned 
#    another way of stating it is that 1-4 is co-linear on query and reference if intervening 1-2 and 3-4 junctions are artifactual
#           =========1==========4=========
#           ||||||||||          ||||||||||
#           ---------1~~~~~~~~~~4---------
#                     ||||||||||
#           ==========2========3==========
# input:
#   name-sorted SITE_SAM on STDIN, ordered from read 5' to 3' end
# output: 
#   name-sorted SITE_SAM on STDOUT
#   where this script sets BLOCK_N

# initialize reporting
our $script = "check_traversal_delta";
our $error  = "$script error";
my ($nReads, $nSvReads, $nJunctions, $nFailedJunctions, $nBlocks) = (0) x 10;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\our $MIN_TRAVERSAL_DELTA, 'MIN_TRAVERSAL_DELTA');

# constants
use constant {
    S_QNAME             => 0, # SITE_SAM fields
    S_FLAG              => 1,
    S_RNAME             => 2,
    S_POS1              => 3, # 1-based
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    DE_TAG              => 6,
    CS_TAG              => 7,
    XF_TAG              => 8,
    XH_TAG              => 9,
    N_REF_BASES         => 10,
    N_READ_BASES        => 11,
    BLOCK_N             => 12,
    SITE_INDEX1_1       => 13,
    SITE_POS1_1         => 14,
    SITE_HAPS_1         => 15,
    SITE_DIST_1         => 16,
    SITE_INDEX1_2       => 17,
    SITE_POS1_2         => 18,
    SITE_HAPS_2         => 19,
    SITE_DIST_2         => 20,
    SEQ_SITE_INDEX1_2   => 21,
    SEQ_SITE_POS1_2     => 22,
    SEQ_SITE_HAPS_2     => 23,
    IS_END_TO_END       => 24,
    READ_HAS_JXN        => 25,
    TARGET_CLASS        => 26,
    S_SEQ               => 27,
    S_QUAL              => 28,
    #-------------
    SPLIT_TO_S_SEQ    => 29,
    #-------------
    _IS_PAIRED      => 1, # SAM FLAG bits
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
    # -------------
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
    # -------------
    FALSE   => 0,
    TRUE    => 1,
};

# parse input SAM
my ($prevRead, @alns);
while(my $aln = <STDIN>){
    my @aln = split("\t", $aln, SPLIT_TO_S_SEQ);
    my $read = $aln[S_QNAME].":".(($aln[S_FLAG] & _IS_PAIRED and $aln[S_FLAG] & _SECOND_IN_PAIR) ? READ2 : READ1);
    if($prevRead and $prevRead ne $read){
        processRead();
        @alns = ();
    }

    push @alns, \@aln;
    $prevRead = $read;
}
processRead();

# print summary information
printCount(commify($nReads),    'nReads',    'input reads');
printCount(commify($nSvReads),  'nSvReads',  'input reads with >1 alignment');
printCount(commify($nJunctions),'nJunctions','junctions in SV reads');
printCount(commify($nFailedJunctions), 'nFailedJunctions', 'junctions that failed traversal delta');
printCount(commify($nBlocks),   'nBlocks',  'output blocks');

# process QNAME alignment groups
sub processRead {
    $nReads++;
    my $blockN = 1;

    # first, 5'-most read alignment always passed as is, i.e., on block 1 (the default)
    print join("\t", @{$alns[0]});

    # only SV-containing reads have additional junction:alignment(s) to process
    if(@alns > 1){
        $nSvReads++;
        $nJunctions += @alns - 1;

        # assess traversal delta for every pair of alignments flanking every possible read sub-path
        # fail all junctions within a failed path
        # once a junction fails, it stays failed even if it passes a different pair of alignments
        my $readLen = length($alns[0][S_SEQ]);
        my @failedTraversalDelta = (0) x @alns;
        foreach my $aln2I(1..$#alns){
            foreach my $aln1I(0..($aln2I - 1)){
                failedTraversalDelta($alns[$aln1I], $alns[$aln2I], $readLen) and
                    @failedTraversalDelta[($aln1I + 1)..$aln2I] = (1) x ($aln2I - $aln1I);
            }
        }

        # print edges with updated BLOCK_N based on failedTraversalDelta
        foreach my $aln2I(1..$#alns){
            $nFailedJunctions += $failedTraversalDelta[$aln2I];
            $failedTraversalDelta[$aln2I] or $blockN++; # thus, increment alignment blockN on every non-failed junction
            $alns[$aln2I][BLOCK_N] = $blockN;
            print join("\t", @{$alns[$aln2I]});
        }
    }
    $nBlocks += $blockN;  
}

# report FALSE if traversal delta is consistent with a true SV, TRUE if failedTraversalDelta
#           =========1==========2========= reference locus flanking the query span
#           ||||||||||          ||||||||||
#           ---------1~~~~~~~~~~2--------- query read, where lower quality ~~ bases are either an unmapped insertion or improperly aligned 
#                     ||||||||||
#           ==========2========1========== a chain of 0 to N false alignments
sub failedTraversalDelta { 
    my ($aln1, $aln2, $readLen) = @_;
    $$aln1[S_RNAME] ne $$aln2[S_RNAME] and return FALSE; # translocation, always a passing traversal delta
    my $strand1 = ($$aln1[S_FLAG] & _REVERSE);
    my $strand2 = ($$aln2[S_FLAG] & _REVERSE);
    $strand1 != $strand2 and return FALSE; # inversion, always a passing traversal delta
    my $refPos1_1 = $strand1 ? $$aln1[S_POS1] : getEnd(@{$aln1}[S_POS1, S_CIGAR]); # i.e., the inner node positions flanking the junction(s)
    my $refPos1_2 = $strand2 ? getEnd(@{$aln2}[S_POS1, S_CIGAR]) : $$aln2[S_POS1]; 
    my $qryPos1_1 = getQueryEnd1($$aln1[S_FLAG], $$aln1[S_CIGAR], $readLen);
    my $qryPos1_2 = getQueryStart0($$aln2[S_FLAG], $$aln2[S_CIGAR]) + 1;
    my $refTraversal = $strand1 ?
        $refPos1_1 - $refPos1_2 :
        $refPos1_2 - $refPos1_1;
    my $queryTraversal = $qryPos1_2 - $qryPos1_1;
    abs($refTraversal - $queryTraversal) < $MIN_TRAVERSAL_DELTA;
}

1;

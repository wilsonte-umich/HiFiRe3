use strict;
use warnings;

# action:
#   calculate the reference vs. read traversal delta between all sets of alignments (across all sets of junctions)
#   use calculated deltas and MIN_TRAVERSAL_DELTA to set BLOCK_N of all alignments
#       any set of alignments where traversalDelta < MIN_TRAVERSAL_DELTA are assigned the same BLOCK_N
#   junctions between alignments with the same BLOCK_N may be considered artifactual during SV analysis
#   block numbering is done before alignment quality truncation when all read alignments are still present
#       unlike other SV metrics, block numbering can depend on alignments distant from the two flanking a junction if nAln > 2
#       we want to ensure the even lower quality alignments have an opportunity to question a read's traversal delta
#   traversal deltas will fail and suppress false junctions when:
#       a large, low quality base string is called as a deletion SV with an equally large insertion
#       a similar base string is falsely aligned to another genomic location when it should have aligned inline
#           =========*==========*========= reference locus 1
#           ||||||||||          ||||||||||
#           ---------*~~~~~~~~~~*--------- query read, where lower quality ~~ bases are either an unaligned insertion or improperly aligned 
#                     ||||||||||
#           ============================== reference locus 2
#       notice that the traversal on reference and query are ~the same from * to *
#       contrast the above with a true deletion with some smaller number of bases inserted at the junction
#           =========*==========*========= reference locus 1
#           ||||||||||          ||||||||||
#           ---------*    ~~    *--------- query read, where ~~ bases are a junction insertion derived from the joining mechanism
#    another way of stating it is that 1-4 is co-linear on query and reference if intervening 1-2 and 3-4 junctions are artifactual
#    such that the first and last segments would not call a SV if the read were aligned end-to-end to reference
#           =========1==========4=========
#           ||||||||||          ||||||||||
#           ---------1~~~~~~~~~~4---------
#                     ||||||||||
#           ==========2========3=========

# report FALSE if traversal delta is consistent with a true SV, TRUE if failedTraversalDelta
#           =========1==========2========= reference locus flanking the query span
#           ||||||||||          ||||||||||
#           ---------1~~~~~~~~~~2--------- query read, where lower quality ~~ bases are either an unmapped insertion or improperly aligned 
#                     ||||||||||
#           ==========2========1========== a chain of 0 to N false alignments

# variables
use vars qw (
    $MIN_TRAVERSAL_DELTA
);

# constants
use constant {
    S_QNAME             => 0, # SITE_SAM fields
    S_FLAG              => 1,
    S_RNAME             => 2,
    S_POS1              => 3,
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    # ... unused fields ...
    # -------------
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
    #-------------
    FALSE   => 0,
    TRUE    => 1,
};

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

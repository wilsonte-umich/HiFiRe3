#!/usr/bin/perl
use strict;
use warnings;

# smith_waterman() performs either Needleman-Wunsch or Smith-Waterman 
# alignment of a query to a reference sequence, depending on the arguments
# expects all upper case IUPAC codes
# returns just one best alignment as an array of query values
#   M operations carry the qry base in the array slot (could be a base mismatch)
#   I operations carry the inserted base prepended to the NEXT reference postion
#   D operations carry "-" in place of the query base that was deleted relative to reference
# in fast mode, query and reference are assumed to be nearly (but not exactly) identical and of similar length
#   such as when comparing duplicate reads of the same sequence span
#   thus, only a limited subset of all possible alignment registers are considered, set by $maxShift
# in forceQryEnd mode, alignment and score must always go to one (but not the other) end of query
#   such as when placing a clip terminus into a candidate region on the other side of an SV junction
#   the end of query to use is set by value of $forceQryEnd, either QRY_START or QRY_END
#   this mode (and only this mode) will fail if there are multiple equally good alignments
# in local mode, best alignments can be local/partial, i.e., they do not need to go to the end of either molecule
# when not forceQryEnd or local, alignments are inclusive to the end of at least one input sequence on each side
# usage: 
#    our ($matchScore, $mismatchPenalty, $gapOpenPenalty, $gapExtensionPenalty, $maxShift) = 
#        (1,           -1.5,             -2.5,            -1,                   3);
#    require "/path/to/smith_waterman.pl";
#    my ($qryOnRef, $score, $startQry, $endQry, $startRef, $endRef) = smith_waterman($qry, $ref, 1); # for fast mode
#    print $score, "\n";
#    print join(" ", @$qryOnRef), "\n";
# options default to those above if not shared with "our"
# with the defaults, max possible score is the length of the shorter input sequence

# load dependencies
use File::Basename;
my $scriptDir = dirname(__FILE__);
require "$scriptDir/IUPAC.pl"; # IUPAC.pl must be in same directory as this script (files split for readability)

# constants for matrix traversal
use constant {
    _DIAG => 0,
    _UP => 1,
    _LEFT => 2,
    QRY_START => 1,
    QRY_END => 2
};

# set scoring metrics to user provided values or defaults
sub setOption {
    my ($score, $default, $negate) = @_;
    defined $$score or $$score = $default;
    $negate and $$score > 0 and $$score = -$$score;  # ensure user-defined penalties are negative
}

# create the lookup table for match/mismatch score for all possible IUPAC code combinations
use vars qw($matchScore $mismatchPenalty $gapOpenPenalty $gapExtensionPenalty $maxShift
            %baseMatches %ryswkmMatches %nMatches
            $suppressEarlyChecks); # used by svWGS to suppress a couple of lines...
sub initializePairScores {
    my %scores;
    foreach my $key(keys %baseMatches)  { $scores{$key} = $matchScore } # e.g. A:A, full match
    foreach my $key(keys %ryswkmMatches){ $scores{$key} = $matchScore / 2 } # e.g. A:R, half match
    foreach my $key(keys %nMatches)     { $scores{$key} = 0 } # e.g. A:N, uninformative base position, neither promoted nor penalized
    $scores{'--'} = 0;
    foreach my $base1(qw (A C G T R Y S W K M N), '-'){
        foreach my $base2(qw (A C G T R Y S W K M N), '-'){
            my $key = "$base1$base2";
            defined $scores{$key} or $scores{$key} = $mismatchPenalty; # everything else defaults to a mismatch, with its penalty
        }
    }
    return %scores;
}
# initialize scoring parameters
setOption(\$matchScore,           1,   0);
setOption(\$mismatchPenalty,     -1.5, 1);
setOption(\$gapOpenPenalty,      -2.501, 1); # 0.001 ajustment gives slight preference to not opening a single-base terminal gap
setOption(\$gapExtensionPenalty, -1,   1);
setOption(\$maxShift,             3,   0);
our %pairedBaseScores = initializePairScores();

# the Smith-Waterman algorithm itself
sub smith_waterman {

    # collect sequence inputs
    my ($qry, $ref, $fast, $forceQryEnd, $local) = @_; # setting $fast to truthy enforces register shift limitations
    $qry or die "smith_waterman error: missing qry sequence\n";
    $ref or die "smith_waterman error: missing ref sequence\n";
    if($qry eq $ref){ # shortcut the process if seqs are identical
       my $n = length($qry);
       return ([split('', $qry)], $matchScore * $n, 0, $n-1, 0, $n-1);  
    }
    my @qry = split('', $qry);
    my @ref = split('', $ref);
    my $nQ = @qry;
    my $nR = @ref;
    if($forceQryEnd and $forceQryEnd == QRY_START){
        @qry = reverse(@qry); # temporarily reverse sequence to allow same code for forceQryEnd QRY_START and QRY_END
        @ref = reverse(@ref); # algorithm is written for QRY_END (i.e. the right side of qry)
    }

    # fill score matrix based on matches and gaps
    my ($bestScore, @matrix, $bestPath, @paths) = (0);
    for(my $refI = 1; $refI <= $nR; $refI++) {    
        my $refBase = $ref[$refI-1];        
        my $isTermI = $refI == $nR;
        my $minQryI = $fast ? $refI-$maxShift : 1;
        my $maxQryI = $fast ? $refI+$maxShift : $nQ;
        for(my $qryI = $minQryI; $qryI <= $maxQryI; $qryI++) { # limit the possible query to reference base matches
            ($qryI < 1 or $qryI > $nQ) and next;
            my $diag_score = ($matrix[$refI-1][$qryI-1][0] || 0) +  # score of prior path
                             $pairedBaseScores{$refBase.$qry[$qryI-1]};
            my $up_score   =  ($matrix[$refI-1][$qryI][0] or 0) + 
                             (($matrix[$refI-1][$qryI][1] or 0) == _UP   ? $gapExtensionPenalty : $gapOpenPenalty); # gap penalties
            my $left_score =  ($matrix[$refI][$qryI-1][0] or 0) + 
                             (($matrix[$refI][$qryI-1][1] or 0) == _LEFT ? $gapExtensionPenalty : $gapOpenPenalty);           
            my ($score, $pointer) = ($diag_score >= $up_score and $diag_score >= $left_score) ? ($diag_score, _DIAG) :
                                    ($up_score >= $diag_score and $up_score >= $left_score)   ? ($up_score,   _UP) :
                                                                                                ($left_score, _LEFT);
            if($local or $forceQryEnd){ 
                $matrix[$refI][$qryI] = $score > 0 ? [$score, $pointer] : []; # allow trimming of non forced query end
                if($local or $qryI == $nQ){ # ensure that all reported alignments go to end of query, unless local
                    if($score > $bestScore){
                        $bestScore = $score; 
                        @paths = ([$refI, $qryI]);   
                    } elsif($score == $bestScore){
                        push @paths, [$refI, $qryI];  # array of equally good paths
                    }
                } 
            } else { # general untrimmed alignment when requiring end-to-end alignment, e.g. in consensus building
                $matrix[$refI][$qryI] = [$score, $pointer]; # best productive path to this residue pair
                if(($isTermI or $qryI == $nQ) and $score > $bestScore){
                    $bestScore = $score;
                    $bestPath = [$refI, $qryI]; # just keep the first encountered path with the best score 
                }                
            }   
        }
    }
    
    # trace backwards to deconvolute best matching path(s) and alignment map(s)
    if($local or $forceQryEnd){ # demand just one best hit in forceQryEnd mode
        !$local and @paths > 1 and return ([('!') x $nQ], 0);
        $bestPath = $paths[0];
    }
    !$bestPath and return ([('!') x $nQ], 0);    
    my ($maxRefI, $maxQryI, @qryOnRef) = @$bestPath;
    my ($refI, $qryI) = ($maxRefI, $maxQryI);
    while (1) {
        my $pointer = $matrix[$refI][$qryI][1];
        defined $pointer or last;  # occurs either just after leftmost unaligned end or at beginning of a sequence 
        if ($pointer == _DIAG) {
            unshift @qryOnRef, $qry[$qryI-1]; # M operation relative to reference (could be a mismatched base)
            $refI--; $qryI--;
        } elsif ($pointer == _LEFT) { # I operation relative to reference, append to NEXT reference base
            $qryOnRef[0] = $qry[$qryI-1].($qryOnRef[0] || ""); # || in case it is the first time through
            $qryI--;
        } else {  # up = D operation relative to reference, pad with a dummy character
            unshift @qryOnRef, '-';
            $refI--;
        }
    }

    # return best alignment
    if($forceQryEnd and $forceQryEnd == QRY_START){
        @qryOnRef = reverse(@qryOnRef); # revert back to original sequence orientation when QUERY_START
        $maxQryI--;
        $maxRefI--;
        $qryI    = $nQ - $qryI - 1;
        $maxQryI = $nQ - $maxQryI - 1;
        $refI    = $nR - $refI - 1;
        $maxRefI = $nR - $maxRefI - 1;
        ($qryI, $maxQryI, $refI, $maxRefI) = ($maxQryI, $qryI, $maxRefI, $refI);
        $maxQryI++;
        $maxRefI++;
    } 
    return (\@qryOnRef, # the output array of all query values 
            $bestScore, # the alignment score             
            $qryI, $maxQryI-1, $refI, $maxRefI-1); # 0-based positions of alignment endpoints on qry and ref         
}

1;

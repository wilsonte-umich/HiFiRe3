use strict;
use warnings;

#========================================================================
# 'exclude.pl' provides functions to suppress alignments/features in bad/excluded genome regions
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
use vars qw($error %chromIndex);
my (@tmpExcludedRegions, @excludedRegions);
my $exclusionsBinSize = 1e6; # split genome into 1 Mb bins for faster searching (but higher memory use)
#========================================================================

#========================================================================
# collect and process the excluded, i.e., reference, features
#------------------------------------------------------------------------
sub initializeExclude {
    my ($excludefile, $nPrependBases) = @_;
    $excludefile or return;
    $nPrependBases or $nPrependBases = 0;
    loadExclusionsBedFile($excludefile, \@tmpExcludedRegions, $nPrependBases);
    collapseExcludedRegions();
    @tmpExcludedRegions = ();
    binExcludedRegions();
}

#========================================================================
# check whether a query feature should be excluded
#------------------------------------------------------------------------
sub isExcludedSpan {  
    my ($chromI, $testStart0, $testEnd1) = @_;
    $chromI or return 1; # exlude non-canonical chromosomes
    $excludedRegions[$chromI] or return 0; # no excluded regions on the chromosome
    $testEnd1 - 1 < $testStart0 and ($testStart0, $testEnd1) = ($testEnd1 - 1, $testStart0 + 1);
    my ($startBinI, $endBinI) = getExcludeCrossedBins($testStart0, $testEnd1); 
    for (my $binI = $startBinI; $binI <= $endBinI; $binI++){ # check all bins crossed by query feature
        $excludedRegions[$chromI][$binI] or next; # no excluded regions in this chromosome bin
        foreach my $refRegion(@{$excludedRegions[$chromI][$binI]}){  
            $$refRegion[0] <= $testEnd1 and $$refRegion[1] >= $testStart0 and return 1; # any overlap with exclude region excludes the test feature            
        }
    }
    return 0;
}
sub isExcludedPosition {
    my ($chromI, $testPos1) = @_;
    $chromI or return 1; # exlude non-canonical chromosomes
    $excludedRegions[$chromI] or return 0; # no excluded regions on the chromosome
    my $binI = int(($testPos1 - 1) / $exclusionsBinSize);
    $excludedRegions[$chromI][$binI] or return 0; # no excluded regions in this chromosome bin
    foreach my $refRegion(@{$excludedRegions[$chromI][$binI]}){
        $$refRegion[0] <= $testPos1 and $$refRegion[1] >= $testPos1 and return 1;
    }
    return 0;
}
#========================================================================

#========================================================================
# worker subs
#------------------------------------------------------------------------
sub loadExclusionsBedFile {
    my ($bedFile, $regions, $nPrependBases) = @_;
    my $inH;
    if($bedFile =~ m/\.gz$/){
        open $inH, "-|", "zcat $bedFile" or die "$error: could not open $bedFile: $!\n";
    } else {
        open $inH, "<", $bedFile         or die "$error: could not open $bedFile: $!\n";
    }
    while (my $line = <$inH>){
        $line =~ m|^\s*#| and next; # ignore comment lines
        chomp $line;
        $line =~ s/\r//g;
        my ($chrom, $start, $end) = split("\t", $line, 4);
        $chrom or next; # ignore blank lines
        my $chromI = $chromIndex{$chrom} or next; # ignore exluded chroms
        parseExcludeInt(\$start); # start and end must be present and integer numbers
        parseExcludeInt(\$end);  
        $start = $start - $nPrependBases; # pad excluded regions at the left end for isExcludedPos
        $start < 0 and $start = 0;
        $$regions[$chromI]{$start}{$end}++;
    }  
    close $inH;
}
sub parseExcludeInt { # validate and uncommify integer values
    my ($int) = @_; 
    defined $$int or die "$error: missing start or end value in exclusions file\n";
    $$int =~ s/,//g;
    $$int =~ m|\D| and die "$error: invalid integer in exclusions file: $int\n";
}
#------------------------------------------------------------------------
sub collapseExcludedRegions { # parse to non-overlapping exclusion regions
    foreach my $chromI(1..$#tmpExcludedRegions){   
        $tmpExcludedRegions[$chromI] or next; 
        my ($regionEnd, $regionStart)= (0); 
        foreach my $start(sort {$a <=> $b} keys %{$tmpExcludedRegions[$chromI]}){     
            if(!$regionEnd){
                $regionStart = $start; 
            } elsif($start > $regionEnd){
                push @{$excludedRegions[$chromI]}, [$regionStart, $regionEnd];
                ($regionStart, $regionEnd) = ($start, 0);
            } 
            foreach my $end(keys %{$tmpExcludedRegions[$chromI]{$start}}){
                $regionEnd >= $end or $regionEnd = $end;
            }
        }
        push @{$excludedRegions[$chromI]}, [$regionStart, $regionEnd];      
    } 
}
#------------------------------------------------------------------------
sub binExcludedRegions {
    my @binned;
    foreach my $chromI(1..$#excludedRegions){  
        $excludedRegions[$chromI] or next; 
        foreach my $refRegion(@{$excludedRegions[$chromI]}){
            my ($startBinI, $endBinI) = getExcludeCrossedBins(@$refRegion);
            for (my $binI = $startBinI; $binI <= $endBinI; $binI++){
                push @{$binned[$chromI][$binI]}, $refRegion;
            }
        }
    } 
    @excludedRegions = @binned;
}
sub getExcludeCrossedBins { # the lowest and highest binIs crossed by a region
    my ($start0, $end1) = @_;
    return (int($start0     / $exclusionsBinSize),
            int(($end1 - 1) / $exclusionsBinSize)); # end bins are converted to 0-referenced start coordinates
}
#========================================================================

1;

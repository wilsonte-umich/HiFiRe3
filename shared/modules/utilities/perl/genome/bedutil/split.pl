#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'split' takes feature BED lines on STDIN and splits any overlapping 
# features into new smaller features corresponding to each unique overlap 
# segment.  Thus, the following four input features:
#    |------------|                    |-----|
#           |-------------|   |----|
# would be split as follows:
#    |------|-----|                    |-----|
#           |-----|-------|   |----|
# or as follows when GROUP_BY=TRUE:
#    |------|-----|-------|   |----|   |-----|
# or as follows when GROUP_BY=TRUE and INCLUDE_GAPS=TRUE:
# |--|------|-----|-------|---|----|---|-----|---------------|
#----------------------------------------------------------------------
# Options are:
#     STRANDED          boolean indicating whether features are strand-specific
#                       if stranded, only features on the same strand are considered overlapping
#                       OPTIONAL [default: FALSE] 
#     INCLUDE_GAPS      boolean indicating whether to include gap features between the split features
#                       input features extending past the end of a chromosome are truncated
#                       OPTIONAL [default: FALSE, only the split features themselves are returned] 
#     CHROM_FILE        file with chromosome information
#                       format = "CHROM\tSIZE", one chromosome per line
#                       REQUIRED if INCLUDE_GAPS=TRUE
#     USE_NONCANONICAL  boolean indicating whether to include non-canonical chromosomes in gap regions
#                       a non-canonical chromosome is anything other than: chr(\d+), chrX, chrY or chrM
#                       OPTIONAL [default: FALSE] 
#     USE_MITO          boolean indicating whether to include the mitochondial chrM in gap regions
#                       OPTIONAL [default: FALSE]  
#     GROUP_BY          boolean indicating whether to group the split output features by coordinate
#                       if not grouped, split features retain the name and score of their parent
#                       OPTIONAL [default: FALSE, duplicate split features are all returned] 
#----------------------------------------------------------------------
# The following options are inherited from group:
#     SCORE_TYPE        how to aggregate scores into a single output score
#                       recognized score types are summarized below
#                       OPTIONAL [default: COUNT] 
#     NAME_JOIN_CHAR    this string is used to join input feature names into the output name
#                       OPTIONAL [default: ",", e.g. "name1,name2,name3"] 
#     FORCE_NAME        use this value as the single output name instead of joining names
#                       OPTIONAL [default: "", i.e. joined names are used] 
#----------------------------------------------------------------------
# The following score types are inherited from group:
#     COUNT             the number of input features that included the split output feature
#     SUM               the sum of scores of the input features that included the split output feature
#     AVERAGE           the average score of the input features that included the split output feature
#     MEDIAN            the median score of the input features that included the split output feature
#     MIN               the minimum score of the input features that included the split output feature
#     MAX               the maximum score of the input features that included the split output feature
#     BOOLEAN           score=1 for any split feature, 0 for any gap between the split features
#---------------------------------------------------------------------- 
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>  
# with no spaces, for example:
#     bedutil split INCLUDE_GAPS CHROM_FILE=/path/to/file
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 60, 54);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    STRANDED
    INCLUDE_GAPS
    CHROM_FILE
    USE_NONCANONICAL
    USE_MITO
    GROUP_BY
    SCORE_TYPE
    NAME_JOIN_CHAR
    FORCE_NAME
);
%booleanOptions = map { $_ => 1} qw (
    STRANDED
    INCLUDE_GAPS
    USE_NONCANONICAL
    USE_MITO
    GROUP_BY
);
%defaultValues = (
    CHROM_FILE => "",
    SCORE_TYPE => "COUNT",
    NAME_JOIN_CHAR => ",",
    FORCE_NAME => ""
);
@requiredOptions = qw (
);
my ($error, $feedback) = parseOptions($scriptName);
#######################################################################


#######################################################################
# collect the chromosome information
#----------------------------------------------------------------------
my ($nChroms, %chromSizes) = (0, 0, 0);
if($ENV{INCLUDE_GAPS}){
    $ENV{CHROM_FILE} or die "$error: missing CHROM_FILE\n";
    open my $chromH, "<", $ENV{CHROM_FILE} or die "$error: could not open $ENV{CHROM_FILE}: $!\n";
    while (my $line = <$chromH>){
        $line =~ m|^\s*#| and next;  # ignore comment lines
        chomp $line;
        $line =~ s/\r//g;
        my ($chrom, $chromSize) = split("\t", $line);
        $chrom or next;  # ignore blank lines
        my $ucChrom = uc($chrom);
        $ENV{USE_NONCANONICAL} or ($ucChrom =~ m/CHR(M|X|Y|\d+)/ or next); 
        $ENV{USE_MITO} or ($ucChrom eq "CHRM" and next);          
        my $error = "$error: invalid chromosome:\n$line\n";  
        parseInt(\$chromSize, $error);        
        $chromSizes{$chrom} = $chromSize;  # chrom sizes needed for last gap, when requested
        $nChroms++; 
    }
    close $chromH;
    $nChroms or die "$error: no chromosomes found in $ENV{CHROM_FILE}\n";
    print STDERR "$feedback: $nChroms chromosomes extracted from $ENV{CHROM_FILE}\n";
}
#######################################################################


#######################################################################
# collect and collapse BED features from STDIN
#----------------------------------------------------------------------
$ENV{IS_CHILD} and return 1;
$ENV{IS_CHILD} = 1;
require "$scriptDir/group.pl";
#----------------------------------------------------------------------
my $gapName = "GAP";
loadBedFile(undef, $error, \my$nFeatures, undef, \my%inFeatures, 1);
$nFeatures or exit;
print STDERR "$feedback: $nFeatures features provided on STDIN\n";
splitBED(\%inFeatures, \my%splitFeatures);
%inFeatures = ();
if($ENV{GROUP_BY}){
    groupBED(\%splitFeatures, \my%groupedFeatures, $gapName);
    %splitFeatures = ();
    printFeaturesHash(*STDOUT, \%groupedFeatures);
} else {
    printFeaturesHash(*STDOUT, \%splitFeatures);
}
#######################################################################


#######################################################################
# split input features
# $features{$chrom}{$strandIndex}{$start}{$end} = [[trailing fields], ...]
#----------------------------------------------------------------------
sub splitBED {  # split overlapping features into multiple features
    my ($inFeatures, $splitFeatures) = @_;
    my @chroms = $ENV{INCLUDE_GAPS} ? keys %chromSizes : keys %$inFeatures;    
    foreach my $chrom(@chroms){
        my @strandIndexes = $ENV{INCLUDE_GAPS} ? ($ENV{STRANDED} ? ("+", "-") : (0)) : keys %{$$inFeatures{$chrom}};
        foreach my $strandIndex(@strandIndexes){
            splitChromStrand($chrom, $strandIndex, $inFeatures, $splitFeatures);
        }
    }
}
sub splitChromStrand {  # split a chromosome on a strand
    my ($chrom, $strandIndex, $inFeatures, $splitFeatures) = @_;
    my $inF = $$inFeatures{$chrom}{$strandIndex};
    unless(scalar(keys %$inF)){
        commitGapFeature($chrom, $strandIndex, 0, $chromSizes{$chrom}, $splitFeatures);
        return;
    } 
    my %splitEnds;  # collect all required end points in the chromosome/strand
    foreach my $start(keys %$inF){
        $start and $splitEnds{$start} = 1;  # a non-zero start implies an end immediately preceding
        foreach my $end(keys %{$$inF{$start}}){
            $splitEnds{$end} = 1;  # ends are taken as is
        }
    }
    my @splitEnds = sort {$a <=> $b} keys %splitEnds;  # ends of _split_ features
    my $i = 0;  # split end counter
    my $maxEnd = 0;  # keep track of the last chromosome position covered by features
    foreach my $start(sort {$a <=> $b} keys %$inF){   # starts of _unsplit_ features   
        commitGapFeature($chrom, $strandIndex, $maxEnd, $start, $splitFeatures);
        while($splitEnds[$i] <= $start){ $i++ }
        my @ends = sort {$a <=> $b} keys %{$$inF{$start}};
        foreach my $end(@ends){
            my $trailings = $$inF{$start}{$end}; # process each group of common-end features
            splitRecursively($chrom, $strandIndex, $start, $end, \@splitEnds, $i, $trailings, $splitFeatures);   
        }
        my $lastEnd = pop @ends;  
        $maxEnd >= $lastEnd or $maxEnd = $lastEnd;
    }
    my $splitStart = pop @splitEnds;  # add gap to end of chromosome
    commitGapFeature($chrom, $strandIndex, $splitStart, $chromSizes{$chrom}, $splitFeatures);
}
sub commitGapFeature {  # if requested, add features for chromosome gaps
    my ($chrom, $strandIndex, $splitStart, $splitEnd, $splitFeatures) = @_;
    $ENV{INCLUDE_GAPS} or return;
    defined $splitEnd or die "$error: unknown chromosome: $chrom\n";    
    $splitEnd > $splitStart or return;
    commitOutFeatures($chrom, $strandIndex, $splitStart, $splitEnd, [[$gapName, 0, $strandIndex?$strandIndex:"+"]], $splitFeatures); 
}
sub commitOutFeatures {  # print an output feature for each input feature crossing the split feature
    my ($chrom, $strandIndex, $splitStart, $splitEnd, $trailings, $splitFeatures) = @_;
    $ENV{INCLUDE_GAPS} and $splitEnd > $chromSizes{$chrom} and $splitEnd = $chromSizes{$chrom};  # truncate features as needed
    foreach my $trailing(@$trailings){
        push @{$$splitFeatures{$chrom}{$strandIndex}{$splitStart}{$splitEnd}}, $trailing;
    }  
}
sub splitRecursively {  # iteratively split the input feature group until split end is at the group end
    no warnings 'recursion';
    my ($chrom, $strandIndex, $splitStart, $end, $splitEnds, $i, $trailings, $splitFeatures) = @_;
    commitOutFeatures($chrom, $strandIndex, $splitStart, $$splitEnds[$i], $trailings, $splitFeatures); 
    if($$splitEnds[$i] < $end){  # features are being split
        $splitStart = $$splitEnds[$i];
        $i++;    
        splitRecursively($chrom, $strandIndex, $splitStart, $end, $splitEnds, $i, $trailings, $splitFeatures);   
    }     
}
#######################################################################

1;

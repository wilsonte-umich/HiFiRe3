#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'group' takes feature BED lines on STDIN and groups them by coordinate.
# Thus, the following six input features:
#    |------|-----|                    |-----|
#           |-----|-------|   |----|
# would be grouped as follows:
#    |------|-----|-------|   |----|   |-----|
#----------------------------------------------------------------------
# Options are:
#     STRANDED          boolean indicating whether features are strand-specific
#                       if stranded, only features on the same strand are grouped
#                       OPTIONAL [default: FALSE] 
#     SCORE_TYPE        how to aggregate scores into a single output score
#                       recognized score types are summarized below
#                       OPTIONAL [default: COUNT] 
#     NAME_JOIN_CHAR    this string is used to join input feature names into the output name
#                       OPTIONAL [default: ",", e.g. "name1,name2,name3"] 
#     FORCE_NAME        use this value as the single output name instead of joining names
#                       OPTIONAL [default: "", i.e. joined names are used] 
#----------------------------------------------------------------------
# Score types are:
#     COUNT             the number of input features that were grouped
#     SUM               the sum of scores of the grouped features
#     AVERAGE           the average score of the grouped features
#     MEDIAN            the median score of the grouped features
#     MIN               the minimum score of the grouped features
#     MAX               the maximum score of the grouped features
#     BOOLEAN           set all scores to 1, i.e. boolean TRUE
#---------------------------------------------------------------------- 
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>  
# with no spaces, for example:
#     bedutil group STRANDED FORCE_NAME=newName
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 38, 32);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    STRANDED
    SCORE_TYPE
    NAME_JOIN_CHAR
    FORCE_NAME
);
%booleanOptions = map { $_ => 1} qw (
    STRANDED
);
%defaultValues = (
    SCORE_TYPE => "COUNT",
    NAME_JOIN_CHAR => ",",
    FORCE_NAME => ""
);
@requiredOptions = qw (
);
my ($error, $feedback) = parseOptions($scriptName);
#----------------------------------------------------------------------
$ENV{IS_CHILD} and return 1;
#######################################################################


#######################################################################
# collect and collapse BED features from STDIN
#----------------------------------------------------------------------
loadBedFile(undef, $error, \my$nFeatures, undef, \my%inFeatures);
$nFeatures or exit;
print STDERR "$feedback: $nFeatures features provided on STDIN\n";
groupBED(\%inFeatures, \my%groupedFeatures);
%inFeatures = ();
printFeaturesHash(*STDOUT, \%groupedFeatures);
#######################################################################


#######################################################################
# group output features
#----------------------------------------------------------------------
sub groupBED {  # group features by chromosome, strand and coordinate
    my ($inFeatures, $groupedFeatures, $ignoreName) = @_;
    foreach my $chrom(keys %$inFeatures){
        foreach my $strandIndex(keys %{$$inFeatures{$chrom}}){
            groupChromStrand($chrom, $strandIndex, $inFeatures, $groupedFeatures, $ignoreName); 
        }
    }
}
sub groupChromStrand {  # group a chromosome on a strand
    my ($chrom, $strandIndex, $inFeatures, $groupedFeatures, $ignoreName) = @_;
    my $inF = $$inFeatures{$chrom}{$strandIndex};
    foreach my $start(keys %$inF){   
        foreach my $end(keys %{$$inF{$start}}){
            my (@names, $sum, $count, @scores);
            foreach my $trailing(@{$$inF{$start}{$end}}){
                my ($name, $score, $strand) = @$trailing;
                unless(!$ignoreName or $name eq $ignoreName){  # do not allow split GAP score=0 to be overwritten
                    push @names, $name;
                    $sum += $score;
                    $count++;
                    push @scores, $score;
                } 
            }  
            my $name = getGroupName(\@names, $ignoreName);
            my $score = getGroupScore($sum, $count, \@scores);
            my $strand = $strandIndex ? $strandIndex : "+";
            push @{$$groupedFeatures{$chrom}{$strandIndex}{$start}{$end}}, [$name, $score, $strand]; 
        }
    }
}
sub getGroupName {  # interpret grouped feature name
    my ($names, $ignoreName) = @_;
    $ENV{FORCE_NAME} and return $ENV{FORCE_NAME};
    return @$names ? join($ENV{NAME_JOIN_CHAR}, sort {$a cmp $b} @$names) : $ignoreName; 
}
sub getGroupScore {  # return the requested group score
    my ($sum, $count, $scores) = @_;
    $count or return 0;
    if($ENV{SCORE_TYPE} eq 'AVERAGE'){
        return $sum / $count;
    } elsif($ENV{SCORE_TYPE} eq 'SUM'){
        return $sum;
    } elsif($ENV{SCORE_TYPE} eq 'COUNT'){
        return $count;
    } elsif($ENV{SCORE_TYPE} eq 'MEDIAN'){
        my @sorted = sort {$a <=> $b} @$scores;
        my $i = int($count/2+0.5);
        return $sorted[$i];
    } elsif($ENV{SCORE_TYPE} eq 'MIN'){
        my @sorted = sort {$a <=> $b} @$scores;
        return $sorted[0];
    } elsif($ENV{SCORE_TYPE} eq 'MAX'){
        my @sorted = sort {$b <=> $a} @$scores;
        return $sorted[0];
    } elsif($ENV{SCORE_TYPE} eq 'BOOLEAN'){
        return 1;
    } else {
        die "$error: unknown group score type: $ENV{SCORE_TYPE}\n";
    }    
}
#######################################################################

1;

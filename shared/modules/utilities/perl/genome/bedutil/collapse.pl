#!/usr/bin/perl
use strict;
use warnings;
use Cwd qw(abs_path);

#######################################################################
# 'collapse' takes feature BED lines on STDIN and collapses them into
# contiguous overlapping BED features, which are printed to STOUT.
# Thus, the following four input features:
#    |------------|                       |-----|
#           |-------------|       |----|
# would be collapsed as follows:
#    |--------------------|       |----|  |-----|
# or as follows with an appropriate PADDING value:
#    |--------------------|       |-------------|
#----------------------------------------------------------------------
# Options are:
#     PADDING      bp added to each end of each feature prior to collapse.
#                  PADDING > 0 will thus merge adjacent but non-overlapping
#                  features into a single contiguous feature that includes
#                  the padded empty gap space between them.
#                  OPTIONAL [default: 0, i.e unpadded]  
#     UNPAD        boolean whether to remove padding from the flanks of collapsed regions
#                  OPTIONAL [default: 0, i.e collapsed features retain outside padding]  
#     STRANDED     boolean indicating whether features are strand-specific
#                  if stranded, only features on the same strand are collapsed together
#                  if unstranded, the strand field is always set to '+'
#                  OPTIONAL [default: FALSE] 
#     COUNT        boolean indicating whether to count the number of collapsed features
#                  the count is appended as the last output column
#                  OPTIONAL [default: FALSE] 
#---------------------------------------------------------------------- 
# Options may be provided as environment variables or on the command line in format:
#     <value_option>=<value> <boolean_option>  
# with no spaces, for example:
#     bedutil collapse PADDING=1000 STRANDED
#######################################################################


#######################################################################
# detect and check options
#----------------------------------------------------------------------
my $scriptFile = abs_path($0);
my ($scriptDir, $scriptName) = $scriptFile =~ m|(.+)/(.+)|;
require "$scriptDir/common.pl";
showHelp($scriptFile, 36, 31);
#----------------------------------------------------------------------
use vars qw(%validOptions %booleanOptions %defaultValues @requiredOptions);
%validOptions = map { $_ => 1} qw (
    PADDING
    UNPAD
    STRANDED
    COUNT
);
%booleanOptions = map { $_ => 1} qw (
    UNPAD
    STRANDED 
    COUNT
);
%defaultValues = (
    PADDING => 0
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
loadBedFile(undef, $error, \my$nFeatures, undef, \my%features);
$nFeatures or exit;
print STDERR "$feedback: $nFeatures features provided on STDIN\n";
collapseBED(\%features, \my%regions);
printRegionsHash(*STDOUT, \%regions);
#######################################################################


#######################################################################
# collapse input 'features' into output 'regions', where:
#    $features{$chrom}{$strandIndex}{$start}{$end} = [[trailing fields], ...]
#    $regions{$chrom}{$strandIndex} = [[$regionStart, $regionEnd, $regionCount, [$name, ...]], ... ]
#----------------------------------------------------------------------
sub collapseBED {  # merge overlapping and immediately adjacent regions into contiguous regions
    my ($features, $regions) = @_;
    foreach my $chrom(keys %$features){
        foreach my $strandIndex(keys %{$$features{$chrom}}){
            collapseChromStrand($chrom, $strandIndex, $features, $regions);
        }
    }
}
sub collapseChromStrand {
    my ($chrom, $strandIndex, $features, $regions) = @_;
    my ($regionStart, $regionCount, @names); 
    my $regionEnd = 0; 
    foreach my $start(sort {$a <=> $b} keys %{$$features{$chrom}{$strandIndex}}){     
        if(!$regionEnd){
            $regionStart = $start; 
        } elsif($start > $regionEnd){
            push @{$$regions{$chrom}{$strandIndex}}, [$regionStart, $regionEnd, $regionCount, [@names]];
            ($regionStart, $regionEnd, $regionCount) = ($start, 0, 0);
            @names = ();
        } 
        foreach my $end(keys %{$$features{$chrom}{$strandIndex}{$start}}){
            $regionEnd >= $end or $regionEnd = $end;
            foreach my $trailing(@{$$features{$chrom}{$strandIndex}{$start}{$end}}){
                $regionCount++;
                push @names, $$trailing[0] ? $$trailing[0] : "$chrom:$start-$end"."[$strandIndex]";
            }
        }
    }
    push @{$$regions{$chrom}{$strandIndex}}, [$regionStart, $regionEnd, $regionCount, [@names]];
}
#######################################################################

1;

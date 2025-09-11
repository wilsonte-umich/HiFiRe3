#!/usr/bin/perl
use strict;
use warnings;

our (%validOptions, %booleanOptions, %defaultValues, @requiredOptions);

#######################################################################
# detect and check options
#----------------------------------------------------------------------
sub showHelp {
    my ($scriptFile, $headLines, $tailLines) = @_;
    if($ARGV[0] and $ARGV[0] eq '--help'){  # provide help if requested
        print "\n", qx/head -n$headLines $scriptFile | tail -n$tailLines | sed 's|^#||'/, "\n";
        exit;
    } 
}
#----------------------------------------------------------------------
sub parseOptions {
    my ($scriptName, $scriptDir) = @_;
    my ($error, $feedback) = ("$scriptName fatal error", $scriptName);
    while (my $arg = shift @ARGV){
        my ($option, $value) = split("=", $arg);
        $option = uc($option);
        $validOptions{$option} or die "$error: unrecognized option: $option\n";
        if($booleanOptions{$option}){
             defined $value or $value = 1;
             my $isFalse = ($value eq "FALSE" or ($value eq "" or $value eq "0"));
             $value = $isFalse ? 0 : 1;
        }
        defined $value or die "$error: missing value for option $option\n";
        $ENV{$option} = $value;  # command line supercedes environment
    }
    foreach my $option(keys %defaultValues){  # environment supercedes default values
        defined $ENV{$option} or $ENV{$option} = $defaultValues{$option};
    }
    foreach my $option(@requiredOptions){
        defined $ENV{$option} or die "$error: a value is required for option $option\n";
    }
    $scriptDir and $ENV{SCRIPT_DIR} = $scriptDir;
    return ($error, $feedback); 
}
#######################################################################


#######################################################################
# load input BED and validate format
#----------------------------------------------------------------------
sub loadBedFile {
    my ($bedFile, $error, $nFeatures, $featuresArray, $featuresHash, $dieOnZeroLength) = @_;
    my $inH = $bedFile ? undef : *STDIN;  # if no file string is provided, collect data from STDIN
    unless($inH){ # if file string has white space, assume it is a stream; otherwise assume it is a file
        my $mode = $bedFile =~ m|\s+| ? "-|" : "<";
        open $inH, $mode, $bedFile or die "$error: could not open $bedFile: $!\n";
    }
    while (my $line = <$inH>){  
        my ($strandIndex, $feature) = parseBedLine($line, $error, $dieOnZeroLength);
        $feature or next;
        $$nFeatures++;        
        $featuresArray and push @$featuresArray, $feature;
        $featuresHash or next;
        my @feature = @$feature;
        my $nFields = @feature;
        my $trailingFields = ($nFields > 3) ? [@feature[3..$#feature]] : [];
        push @{$$featuresHash{$feature[0]}{$strandIndex}{$feature[1]}{$feature[2]}}, $trailingFields;
    }  
    close $inH;
}
#----------------------------------------------------------------------
sub parseBedLine {
    my ($line, $error, $dieOnZeroLength) = @_;
    $line =~ m|^\s*#| and return;  # ignore comment lines
    chomp $line;
    $line =~ s/\r//g;
    my @feature = split("\t", $line);
    $feature[0] or return;  # ignore blank lines
    $error = "$error: invalid BED feature:\n$line\n";    
    parseInt(\$feature[1], $error);  # start and end must be present and integer numbers
    parseInt(\$feature[2], $error);
    $dieOnZeroLength and ($feature[2] - $feature[1] > 0 or
                          die "zero or negative length in line:\n$line\n");
    checkStrand($feature[5], $error);
    padFeature(\@feature);
    my $strandIndex = $ENV{STRANDED} ? $feature[5] : 0;
    return ($strandIndex, \@feature);
}
#----------------------------------------------------------------------
sub parseInt {  # validate and uncommify integer values
    my ($int, $error) = @_;  # int passed as reference
    $error or $error = "";
    defined $$int or die $error."missing an expected integer\n";
    my $unparsed = $$int;
    $$int =~ s/,//g;
    $$int =~ m|\D| and die $error."invalid integer: $unparsed\n";
}
sub checkStrand {
    my ($strand, $error) = @_;
    $ENV{STRANDED} or return;
    $error or $error = "";
    defined $strand or die $error."BED strand expected when STRANDED is true\n";
    ($strand eq '+' or $strand eq '-') or die $error."invalid strand: $strand\n";
}
#######################################################################


#######################################################################
# handle feature padding
#----------------------------------------------------------------------
sub padFeaturesHash {
    my ($features) = @_;
    $ENV{PADDING} or return $features;
    my %padded;
    foreach my $chrom(keys %$features){   
        foreach my $strandIndex(keys %{$$features{$chrom}}){
            foreach my $start(keys %{$$features{$chrom}{$strandIndex}}){ 
                my $paddedStart = $start - $ENV{PADDING};
                foreach my $end(keys %{$$features{$chrom}{$strandIndex}{$start}}){ 
                    my $paddedEnd = $end + $ENV{PADDING};
                    foreach my $trailing(@{$$features{$chrom}{$strandIndex}{$start}{$end}}){ 
                        push @{$padded{$chrom}{$strandIndex}{$paddedStart}{$paddedEnd}}, $trailing;
                    } 
                }  
            }
        }
    }  
    return \%padded;
}
sub padFeature {
    my ($feature) = @_;
    $ENV{PADDING} or return;
    $$feature[1] -= $ENV{PADDING};
    $$feature[2] += $ENV{PADDING};
}
#----------------------------------------------------------------------
sub unpadFeature {
    my ($feature, $force) = @_;
    $ENV{PADDING} or return;
    $force or ($ENV{UNPAD} or return); 
    $$feature[1] += $ENV{PADDING};
    $$feature[2] -= $ENV{PADDING};
}
#######################################################################


#######################################################################
# handle feature offsetting
#----------------------------------------------------------------------
sub offsetFeaturesHash {
    my ($features) = @_;
    $ENV{OFFSET_CONTROL} or return $features;
    my %offset;
    foreach my $chrom(keys %$features){   
        foreach my $strandIndex(keys %{$$features{$chrom}}){
            foreach my $start(keys %{$$features{$chrom}{$strandIndex}}){ 
                my $offsetStart = $start + $ENV{OFFSET_CONTROL};
                foreach my $end(keys %{$$features{$chrom}{$strandIndex}{$start}}){ 
                    my $offsetEnd = $end + $ENV{OFFSET_CONTROL};
                    foreach my $trailing(@{$$features{$chrom}{$strandIndex}{$start}{$end}}){ 
                        push @{$offset{$chrom}{$strandIndex}{$offsetStart}{$offsetEnd}}, $trailing;
                    } 
                }  
            }
        }
    }  
    return \%offset;
}
#######################################################################


#######################################################################
# print BED results to STDOUT
#----------------------------------------------------------------------
sub printRegionsHash {
    my ($outH, $regions, $forceName, @append) = @_; 
    foreach my $chrom(keys %$regions){
        foreach my $strandIndex(keys %{$$regions{$chrom}}){
            foreach my $region(@{$$regions{$chrom}{$strandIndex}}){      
                my ($regionStart, $regionEnd, $regionCount, $names) = @$region; 
                my $name = $forceName ? $forceName : join("::", @$names);   
                my $strand = $ENV{STRANDED} ? $strandIndex : '+';      
                my @feature = ($chrom, $regionStart, $regionEnd, $name, 0, $strand);
                $ENV{COUNT} and push @feature, $regionCount;
                @append and push @feature, @append; 
                unpadFeature(\@feature);
                print $outH join("\t", @feature), "\n";
            }
        } 
    }
}
sub printFeaturesHash {
    my ($outH, $features, @append) = @_;
    foreach my $chrom(keys %$features){   
        foreach my $strandIndex(keys %{$$features{$chrom}}){
            foreach my $start(keys %{$$features{$chrom}{$strandIndex}}){ 
                foreach my $end(keys %{$$features{$chrom}{$strandIndex}{$start}}){ 
                    foreach my $trailing(@{$$features{$chrom}{$strandIndex}{$start}{$end}}){ 
                        print $outH join("\t", $chrom, $start, $end, @$trailing, @append), "\n";
                    } 
                }  
            }
        }
    }  
}
sub printFeaturesArray {
    my ($outH, $features, @append) = @_;
    foreach my $feature(@$features){  
        print $outH join("\t", @$feature, @append), "\n";
    }
}
#######################################################################


#######################################################################
# convert feature data types
#----------------------------------------------------------------------
sub featureHashToRegionHash {
    my ($chrom, $strandIndex, $features, $regions) = @_;
    foreach my $start(keys %{$$features{$chrom}{$strandIndex}}){ 
        foreach my $end(keys %{$$features{$chrom}{$strandIndex}{$start}}){ 
            foreach my $trailing(@{$$features{$chrom}{$strandIndex}{$start}{$end}}){ 
                push @{$$regions{$chrom}{$strandIndex}}, [$start, $end, @$trailing];            
            } 
        }  
    }
}
#######################################################################

1;

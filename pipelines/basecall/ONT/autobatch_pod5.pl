use strict;
use warnings;

# action:
#   construct a set of POD5 file batch to use POD_BUFFER for:
#       maximal efficiency
#       no buffer overruns
# input:
#   `ls -l`` or `tar -tv`` of the set of input POD5 files
# ouput:
#   $POD5_BUFFER_DIR/batches.txt
# return value:
#   the number of pod5 batches found in output file

# get the input column with the file size (different for tar and ls) and output file
my ($sizeField0, $batchIndexFile) = @ARGV;

# parse the maximum POD5 buffer size; must be in G units
my $pod5BufferSize = $ENV{POD5_BUFFER_SIZE} or die "missing value for POD5_BUFFER_SIZE\n";
$pod5BufferSize =~ m/\d+(\D+)$/ and (uc($1) eq "G" or die "POD5_BUFFER_SIZE must be in G units\n");
$pod5BufferSize =~ s/\D//g;
$pod5BufferSize *= 1e9;

# parse the minimum POD5 file size; must be in M units
my $minPod5Size = $ENV{MIN_POD5_SIZE} || "0M";
$minPod5Size =~ m/\d+(\D+)$/ and (uc($1) eq "M" or die "MIN_POD5_SIZE must be in M units\n");
$minPod5Size =~ s/\D//g;
$minPod5Size *= 1e6;

# working parameters
my $maxBatchSize = $pod5BufferSize * 0.8 / 2; # thus, each batch gets up to ~40% of POD5_BUFFER_SIZE
my (%pod5Files, %batched, $nPod5Files);

# collect the files names and sizes
while (my $pod5 = <STDIN>){
    chomp $pod5;
    my @f = split(/\s+/, $pod5);
    if($f[$sizeField0] >= $minPod5Size){
        $pod5Files{$f[$#f]} = $f[$sizeField0];
        $nPod5Files++;
    }
}

# sort the files by size
my @pod5FilesHighLow = sort { $pod5Files{$b} <=> $pod5Files{$a} } keys %pod5Files;
my @pod5FilesLowHigh = reverse @pod5FilesHighLow;

# open the output file handle
open my $outH, ">", $batchIndexFile or die "could not open $batchIndexFile: $!\n";

# run the logic to create file batches
# first batch often take just the single biggest file
# second batch works from the smallest files to accumulate just enough to not exceed buffer
# the process repeats
# alternating between biggest and smallest files provides granularity to ensure
# that no sequential pairing of two batch grows excessively large
# it also ensures that there are no tiny batches
my ($i, $j, $batchN, $prevBatchSize) = (0, 0, 1, 0);
while(1){
    getBatch(\$i, \@pod5FilesHighLow) or last;
    getBatch(\$j, \@pod5FilesLowHigh) or last;
}
sub getBatch {
    my ($k, $pod5s) = @_;
    arePendingFiles() or return;
    my $batchSize = 0;
    while($batchSize < $maxBatchSize and 
          $batchSize + $prevBatchSize < 2 * $maxBatchSize and 
          arePendingFiles()){
        my $pod5 = $$pod5s[$$k];
        my $size = $pod5Files{$pod5};
        $batched{$pod5} or print $outH join("\t", $batchN, (int($size / 1e9 * 100 + 0.5) / 100)."G", $pod5), "\n";
        $batched{$pod5}++;
        $batchSize += $size;
        $$k++;
    }
    $prevBatchSize = $batchSize;
    $batchN++;
    return 1;
}
sub arePendingFiles {
    scalar(keys %batched) < $nPod5Files
}

# return the number of file batches
close $outH;
print $batchN - 1, "\n";

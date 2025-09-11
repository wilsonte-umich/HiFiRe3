use strict;
use warnings;

# change the fastp merged tag to tags carried in QNAME for passing through minimap2 as fastq
# count reads and bases to assess what fastp filtering and merging did to prepare_fastq.pl output
# final QNAME output tag format at this stage is:
#   QNAME:channel:trim5:trim3:mergeLevel:nRead1:nRead2
# where this script sets:
#     mergeLevel, 
#         2 for fastp-merged reads (highest priority on sorting)
#         0 for single, unmerged, or orphaned reads
#             mergeLevel may later be adjusted to 1 if subjected to downstream alignment-guided merging
#     nRead1, the number of read1 bases present in a final merged read (0 for single, unmerged or orphan reads)
#     nRead2, the number of read2 bases present in a final merged read

# initialize reporting
our $action = "adjust_merge_tags";
my ($nInputReads, $nMergedReads, $nUnmergedReads, $nUnmergedEvents, $nPairedEvents) = (0) x 10;
my @baseCounts = (0, 0, 0);

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# constants
use constant {
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
    MERGED  => 3,
};

# run the interleaved pairs
my $nullMergeTags = join(":", 0, 0, 0);
my $lineN = 0;
my ($prevQName, $readN) = ("");
while(my $line = <STDIN>){
    if($lineN % 4 == 0){
        $nInputReads++;

        # split drops trailing whitespace when splitting on whitespace, i.e., chomp occurs implicitly
        my @f = split(" ", $line); 

        # any prior platform-specific trailing data on QNAME line was stripped by align/do/prepare_fastq.pl
        # thus, anything present now was added by fastp
        # QNAME merged_150_15 means that 150bp are from read1, and 15bp are from read2
        if($f[1]){ 
            $nMergedReads++;
            my ($nRead1, $nRead2) = ($f[1] =~ m/merged_(\d+)_(\d+)/); 
            $line = join(":", 
                $f[0],
                2,
                $nRead1,
                $nRead2
            )."\n";

        # QNAME only signifies a single or unmerged read
        } else {
            $nUnmergedReads++;
            $line = join(":", 
                $f[0], # incoming QNAME:channel:trim5:trim3
                $nullMergeTags # null mergeLevel:nRead1:nRead2
            )."\n";
            $prevQName eq $line or $nUnmergedEvents++;
        }
        $readN = $prevQName eq $line ? READ2 : READ1;
        $readN == READ2 and $nPairedEvents++;
        $prevQName = $line;
    } elsif($lineN % 4 == 1){
        $baseCounts[$readN] += length($line) - 1;
    }
    print $line;
    $lineN++;
}

# print summary information
my $nInputEvents    = $nMergedReads + $nUnmergedEvents;
my $nOrphanedReads  = $nUnmergedReads - 2 * $nPairedEvents;
my $percentMerged   = $nMergedReads   / $nInputEvents * 100;
my $percentOrphaned = $nOrphanedReads / $nInputEvents * 100;
printCount(commify($nInputReads),       'nInputReads',       'input reads');
printCount(commify($nInputEvents),      'nInputEvents',      'input events');
printCount(commify($nMergedReads),      'nMergedReads',      'input merged reads == events with merged read pairs');
printCount(commify($nUnmergedReads),    'nUnmergedReads',    'input unmerged reads (single, paired, or orphaned reads)');
printCount(commify($nUnmergedEvents),   'nUnmergedEvents',   'input events with unmerged reads (single reads, unmerged read pairs, or orphaned reads)');
printCount(commify($nPairedEvents),     'nPairedEvents',     'input events with paired reads == read2 count');
printCount(commify($nOrphanedReads),    'nOrphanedReads',    'input read1 without a read2 partner = events that lost (or never had) a read partner');
printCount(roundCount2($percentMerged), 'percentMerged',     'percent of reads pairs merged');
printCount(roundCount2($percentOrphaned),'percentOrphaned',  'percent of reads pairs orphaned');
printCount(commify($baseCounts[READ1]), 'baseCounts[READ1]', 'read1 bases in input events (single, first paired, merged, or orphaned)');
printCount(commify($baseCounts[READ2]), 'baseCounts[READ2]', 'read2 bases in input events (second paired)');
printCount(commify($baseCounts[READ1] + $baseCounts[READ2]), 'baseCounts[EVENT]', 'read1 + read2 bases in input events (nearly all unique)');


# from https://github.com/OpenGene/fastp?tab=readme-ov-file#merge-paired-end-reads

# In the output file, a tag like merged_xxx_yyy will be added to each read name 
# to indicate that how many base pairs are from read1 and from read2, respectively. 
# For example,  
#     @NB551106:9:H5Y5GBGX2:1:22306:18653:13119 1:N:0:GATCAG merged_150_15 
# means that 
#     150bp are from read1, and 15bp are from read2
# fastp prefers the bases in read1 since they usually have higher quality than read2.

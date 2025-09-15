use strict;
use warnings;

# action:
#     prepare FASTQ for alignment from unaligned BAM files
#     propagate all tags through unaligned bam through FASTQ to aligned BAM
#     enforce minimum read length of 250bp
#     tabulate observed insert size distributions, i.e., read lengths, before any downstream filtering or projection
# expects:
#     source $MODULES_DIR/align/set_read_file_vars.sh
# output: 
#     FASTQ on stdout with tags as QNAME extensions

# initialize reporting
our $action = "prepare_fastq";
my ($nInputEvents, $nOutputEvents, $baseCounts) = (0) x 10;
my %insertSizeCounts;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\my $READ_FILES,        'READ_FILES');
fillEnvVar(\my $UNFILTERED_INSERT_SIZES_FILE, 'UNFILTERED_INSERT_SIZES_FILE');
$READ_FILES =~ s/[\n|\r]/ /g;

# constants
use constant {
    _QNAME  => 0,
    _TAGS   => 1,
    _SEQ    => 2,
    _QUAL   => 3,
    _PASSED_SIZE_FILTER => 4,
    #-------------
    MIN_INSERT_SIZE => 250,
    SIZE_PLOT_BIN_SIZE => 250
};

# open and pass the file input handles
foreach my $ubamFile(split(" ", $READ_FILES)){
    open my $inH, "-|", "samtools view $ubamFile | samtools fastq -T -" #  ch,MM,ML
        or throwError("could not open $ubamFile: $!");
    runSingleReads($inH);
    close $inH;
}

# run single reads
# no need to multi-thread since the downstream alignment process is rate limiting
sub runSingleReads {
    my ($inH) = @_;
    while (my $read = getRead($inH)){
        $nInputEvents++; 
        my $len = length($$read[_SEQ]);
        $insertSizeCounts{int($len / SIZE_PLOT_BIN_SIZE + 0.5) * SIZE_PLOT_BIN_SIZE}++;
        if($$read[_PASSED_SIZE_FILTER]){
            $nOutputEvents++;
            my $nameLine = $$read[_TAGS] ? "$$read[_QNAME] $$read[_TAGS]" : $$read[_QNAME];
            print join("\n", $nameLine, $$read[_SEQ], '+', $$read[_QUAL]), "\n";
            $baseCounts += $len;
        }
    }
}

# parse a FASTQ set of 4 lines
sub getRead {
    my ($inH) = @_;
    my $name = <$inH>;
    $name or return; # EOF
    chomp $name;
    my ($qName, $tags) = split(/\s/, $name, 2);
    my $seq     = <$inH>;
    my $discard = <$inH>;
    my $qual    = <$inH>; 
    chomp $seq;
    chomp $qual;
    return [
        $qName, 
        $tags, 
        $seq, 
        $qual, 
        length($seq) >= MIN_INSERT_SIZE
    ];
}

# print summary information
printCount(commify($nInputEvents),  'nInputEvents',  'input events, i.e., reads');
printCount(commify($nOutputEvents), 'nOutputEvents', "output events after size filtering >= " . MIN_INSERT_SIZE);
printCount(commify($baseCounts),    'baseCounts',    'read bases in output events');

# print insert size counts
open my $outH, '>', "$UNFILTERED_INSERT_SIZES_FILE" or throwError("could not open insert sizes file: $!");
my @insertSizes = sort {$a <=> $b} keys %insertSizeCounts;
print $outH join("\n", map { 
    join("\t", 
        "unfiltered",
        $_, 
        $insertSizeCounts{$_}
    )
} @insertSizes), "\n";
close $outH;

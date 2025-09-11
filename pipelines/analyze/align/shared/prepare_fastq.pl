use strict;
use warnings;
use File::Basename;
use File::Copy;

# action:
#     prepare interleaved FASTQ from different input types
#     trim fixed-length short-read platform reads to $FIXED_READ_LENGTH
#        typically, this is done to trim the low-quality 3' "extra base"
#     append dummy channel:trim5:trim3 values to QNAME if not already specified upstream for the platform
# expects:
#     source $MODULES_DIR/align/set_read_file_vars.sh
#     source ${MODULES_DIR}/agFree/set_agFree_vars.sh
#     $FIXED_READ_LENGTH when $IS_FIXED_READ_LENGTH==TRUE
# output: 
#     interleaved FASTQ on stdout

# initialize reporting
our $action = "prepare_fastq";
my ($nInputReads, $nInputEvents, $nOutputEvents) = (0) x 10;
my @baseCounts = (0, 0, 0);

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\my $TMP_DIR_WRK_SMALL,  'TMP_DIR_WRK_SMALL');
fillEnvVar(\my $SEQUENCING_PLATFORM,'SEQUENCING_PLATFORM');
fillEnvVar(\my $READ_PAIR_TYPE,     'READ_PAIR_TYPE');
fillEnvVar(\my $IS_FIXED_READ_LENGTH,'IS_FIXED_READ_LENGTH');
fillEnvVar(\my $FIXED_READ_LENGTH,  'FIXED_READ_LENGTH');
fillEnvVar(\my $READ_FILE_TYPE,     'READ_FILE_TYPE');
fillEnvVar(\my $READ_1_FILES,       'READ_1_FILES');
fillEnvVar(\my $READ_2_FILES,       'READ_2_FILES');
fillEnvVar(\my $MIN_INSERT_SIZE,    'MIN_INSERT_SIZE');
fillEnvVar(\my $DROP_PLATFORM_NAMES,'DROP_PLATFORM_NAMES');
fillEnvVar(\my $EXTERNAL_ONT_BASECALLING,'EXTERNAL_ONT_BASECALLING');
my $isONT = ($SEQUENCING_PLATFORM eq "ONT"); # at present, the only platform subject to read trimming prior to running align.sh
my $addNullTrim = (!$isONT or $EXTERNAL_ONT_BASECALLING);
my $nullTrim = join(":", 0, 0, 0);
$READ_1_FILES =~ s/[\n|\r]/ /g;
$READ_1_FILES =~ s/[\n|\r]/ /g;

# constants
use constant {
    _QNAME  => 0,
    _SEQ    => 1,
    _QUAL   => 2,
    _PASSED_SIZE_FILTER => 3,
    #-------------
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
};

# open and pass the file input handles
if($READ_FILE_TYPE eq "fastq.gz"){
    if($READ_PAIR_TYPE eq "paired"){
        open my $inH1, "-|", "zcat $READ_1_FILES" or throwError("could not open $READ_1_FILES: $!");
        open my $inH2, "-|", "zcat $READ_2_FILES" or throwError("could not open $READ_2_FILES: $!");
        runReadPairs($inH1, $inH2);
        close $inH1;
        close $inH2; 
    } else {
        open my $inH, "-|", "zcat $READ_1_FILES" or throwError("could not open $READ_1_FILES: $!");
        runSingleReads($inH);
        close $inH;
    }
} elsif($READ_FILE_TYPE eq "fastq"){
    if($READ_PAIR_TYPE eq "paired"){
        open my $inH1, "-|", "cat $READ_1_FILES" or throwError("could not open $READ_1_FILES: $!");
        open my $inH2, "-|", "cat $READ_2_FILES" or throwError("could not open $READ_2_FILES: $!");
        runReadPairs($inH1, $inH2);
        close $inH1;
        close $inH2; 
    } else {
        open my $inH, "-|", "cat $READ_1_FILES" or throwError("could not open $READ_1_FILES: $!");
        runSingleReads($inH);
        close $inH;
    }
} elsif($READ_FILE_TYPE eq "unaligned.bam"){
    if($READ_PAIR_TYPE eq "paired"){
        die "$action error: paired reads not yet supported for unaligned.bam input\n";
    } else {
        foreach my $ubamFile(split(" ", $READ_1_FILES)){ # might contain single or interleaved paired reads
            open my $inH, "-|", "samtools view $ubamFile | samtools fastq -" or throwError("could not open $ubamFile: $!");
            runSingleReads($inH);
            close $inH;
        }
    }
} elsif($READ_FILE_TYPE eq "sra"){
    foreach my $sraFile(split(" ", $READ_1_FILES)){ # might contain single or interleaved paired reads
        # when operating in a stream, fastq-dump is preferred
        my $sraName = basename($sraFile);
        my $tmpFile = "$TMP_DIR_WRK_SMALL/$sraName";
        unlink $tmpFile; # in case some partial file pre-exists
        copy($sraFile, $tmpFile) or die "SRA tmp copy failed: $sraName: $!"; # pre-copy SRA file to SSD to speed up fastq-dump
        open my $inH, "-|", "fastq-dump --stdout --split-files $tmpFile" or throwError("could not open $tmpFile: $!");
        if($READ_PAIR_TYPE eq "paired"){
            runReadPairs($inH, $inH);
        } else {
            runSingleReads($inH);
        }
        close $inH;
        unlink $tmpFile;
    }
} else {
    die "unrecognized READ_FILE_TYPE: $READ_FILE_TYPE\n"; 
}

# run single reads
# no need to multi-thread since the downstream alignment process is rate limiting
sub runSingleReads {
    my ($inH) = @_;
    while (my $read = getRead($inH)){
        $nInputEvents++;
        if($$read[_PASSED_SIZE_FILTER]){
            $nOutputEvents++;
            $addNullTrim and $$read[_QNAME] = join(":", $$read[_QNAME], $nullTrim);
            $DROP_PLATFORM_NAMES and $$read[_QNAME] =~ s/.+(:\d+:\d+:\d+)$/\@$nOutputEvents$1/;
            print join("\n", $$read[_QNAME], $$read[_SEQ], '+', $$read[_QUAL]), "\n"; 
            $baseCounts[READ1] += length($$read[_SEQ]);
        }
    }
}

# run paired reads
sub runReadPairs {
    my ($inH1, $inH2) = @_;
    while (my $read1 = getRead($inH1)){
        my $read2 = getRead($inH2);
        $nInputEvents++;
        if($$read1[_PASSED_SIZE_FILTER] and $$read2[_PASSED_SIZE_FILTER]){
            $nOutputEvents++;
            $addNullTrim and $$read1[_QNAME] = join(":", $$read1[_QNAME], $nullTrim);
            $DROP_PLATFORM_NAMES and $$read1[_QNAME] =~ s/.+(:\d+:\d+:\d+)$/\@$nOutputEvents$1/;
            print join("\n", $$read1[_QNAME], $$read1[_SEQ], '+', $$read1[_QUAL]), "\n"; # print interleaved read pairs
            print join("\n", $$read1[_QNAME], $$read2[_SEQ], '+', $$read2[_QUAL]), "\n";
            $baseCounts[READ1] += length($$read1[_SEQ]);
            $baseCounts[READ2] += length($$read2[_SEQ]);
        }
    }
}

# parse a FASTQ set of 4 lines
sub getRead {
    my ($inH) = @_;
    my $name = <$inH>;
    $name or return; # EOF
    $nInputReads++;
    ($name) = split(/\s/, $name); # discard any platform-specific non-name content on FASTQ name line
    my $seq     = <$inH>;
    my $discard = <$inH>;
    my $qual    = <$inH>; 
    chomp $name;
    chomp $seq;
    chomp $qual;
    if($IS_FIXED_READ_LENGTH){
        $seq  = substr($seq,  0, $FIXED_READ_LENGTH);
        $qual = substr($qual, 0, $FIXED_READ_LENGTH);
    }
    return [
        $name, 
        $seq, 
        $qual, 
        length($seq) >= $MIN_INSERT_SIZE
    ];
}

# print summary information
printCount(commify($nInputReads),       'nInputReads',       'input reads');
printCount(commify($nInputEvents),      'nInputEvents',      'input events (e.g., read pairs)');
printCount(commify($nOutputEvents),     'nOutputEvents',     'output events after filtering');
printCount(commify($baseCounts[READ1]), 'baseCounts[READ1]', 'read1 bases in output events');
printCount(commify($baseCounts[READ2]), 'baseCounts[READ2]', 'read2 bases in output events');
printCount(commify($baseCounts[READ1] + $baseCounts[READ2]), 'baseCounts[EVENT]', 'read1 + read2 bases in output events (not all unique given pair overlaps)');

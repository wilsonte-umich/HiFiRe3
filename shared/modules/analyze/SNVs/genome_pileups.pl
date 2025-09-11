use strict;
use warnings;

# generic functions to support `genotype SNVs` and `analyze SNVs`

# constants
use constant {
    ALN_STRAND_POS1     => 0, # strictly speaking, these are alignment types ...
    ALN_VA_TAG          => 1,
    ALN_QS_TAG          => 2,
    ALN_READ_OFFET_5    => 3,
    ALN_X_READ1         => 4,
    ALN_X_READ2         => 5,
    ALN_N_OBSERVED      => 6, # ... with observation counts
    #-------------    
    NO_COVERAGE => "!q0Q0",  # used for no-coverage runs in initial stranded pileups (no-coverage runs omitted from final unstranded pileup)
};

# variables
use vars qw(
    $error $TMP_PILEUP_DIR $SNV_GENOME_PILEUP_PREFIX $SNV_SUMMARY_TABLE_PREFIX
    $chrom $paddedChromIndex1 $chromSize $strandIndex0
    $pileupStartPos1 $lastProcessedI0 @pileup
    $pileupH $snvH
);
my ($childN, $snvAlnsStreamFn);

# create the genome pileup for each strand
sub runStrandedGenomePileup {
    ($childN, $strandIndex0, $snvAlnsStreamFn) = @_;
    my $message = "thread $childN\trunStrandedGenomePileup $chrom $strandIndex0";
    printTimestampedMessage($message); 

    # open file handles
    my $snvAlnsStream = &$snvAlnsStreamFn();
    open my $snvAlnsH, "-|", $snvAlnsStream
        or die "$error: could not open: SNV alns stream: $!\n";
    open $pileupH, ">", "$TMP_PILEUP_DIR/genome_pileup.$paddedChromIndex1.$strandIndex0.txt" # cannot be compressed since need tac below
        or die "$error: could not open: output file: $!\n";

    # run the indexed and sorted alignment groups one at a time
    # alignment groups come in 5' to 3' order on the reference $strandIndex0
    ($pileupStartPos1, $lastProcessedI0, @pileup) = (1, -1);
    while(my $line = <$snvAlnsH>){
        chomp $line;
        my @aln = split("\t", $line);
        $aln[ALN_STRAND_POS1] > $pileupStartPos1 and finishBases_genome($aln[ALN_STRAND_POS1] - 1);
        incrementPileups_genome(\@aln);
    }
    finishBases_genome($pileupStartPos1 + @pileup - 1);
    print $pileupH join("\t", scalar(@pileup), $pileup[0]), "\n";
    my $firstUnfinishedPos1 = $pileupStartPos1 + @pileup;
    $firstUnfinishedPos1 <= $chromSize and print $pileupH join("\t", $chromSize - $firstUnfinishedPos1 + 1, NO_COVERAGE), "\n";

    # finish up 
    printTimestampedMessage($message." done"); 
    close $snvAlnsH;
    close $pileupH;
}

# merge the two strand pileups into a single pileup
sub mergeStrandedGenomePileups {
    ($childN) = @_;
    my $message = "thread $childN\tmergeStrandedGenomePileups $chrom";
    printTimestampedMessage($message); 

    # open file handles
    open my $inH_top, "-|", "cat $TMP_PILEUP_DIR/genome_pileup.$paddedChromIndex1.0.txt" 
        or die "$error: could not open: $!\n";
    open my $inH_bot, "-|", "tac $TMP_PILEUP_DIR/genome_pileup.$paddedChromIndex1.1.txt" # tac, not cat, to work in top-strand order
        or die "$error: could not open: $!\n";
    open $pileupH, "|-", "gzip -c > $SNV_GENOME_PILEUP_PREFIX.$paddedChromIndex1.bed.gz"
        or die "$error: could not open: $!\n";
    open $snvH, "|-", "sort -k2,2n -k3,3n -S 1G | gzip -c > $SNV_SUMMARY_TABLE_PREFIX.$paddedChromIndex1.bed.gz"
        or die "$error: could not open: $!\n";

    # merge the two strand pileups, splitting their runs against each other into smaller runs
    my $refPos1 = 1;
    my ($nPosInRun_top, $runSpec_top, $startPos1_top, $endPos1_top) = getStrandPileup($inH_top, $refPos1);
    my ($nPosInRun_bot, $runSpec_bot, $startPos1_bot, $endPos1_bot) = getStrandPileup($inH_bot, $refPos1);
    NEXT_MERGE_CHUNK: if($endPos1_top and $endPos1_bot){
        if($endPos1_top == $endPos1_bot){
            $refPos1 = mergePileupRuns($refPos1, $nPosInRun_top, $runSpec_top, $runSpec_bot);
            ($nPosInRun_top, $runSpec_top, $startPos1_top, $endPos1_top) = getStrandPileup($inH_top, $refPos1);
            ($nPosInRun_bot, $runSpec_bot, $startPos1_bot, $endPos1_bot) = getStrandPileup($inH_bot, $refPos1);
        } elsif($endPos1_top < $endPos1_bot){
            $refPos1 = mergePileupRuns($refPos1, $nPosInRun_top, $runSpec_top, $runSpec_bot);
            $nPosInRun_bot -= $nPosInRun_top;
            ($nPosInRun_top, $runSpec_top, $startPos1_top, $endPos1_top) = getStrandPileup($inH_top, $refPos1);
        } else {
            $refPos1 = mergePileupRuns($refPos1, $nPosInRun_bot, $runSpec_top, $runSpec_bot);
            $nPosInRun_top -= $nPosInRun_bot;
            ($nPosInRun_bot, $runSpec_bot, $startPos1_bot, $endPos1_bot) = getStrandPileup($inH_bot, $refPos1);
        }
        goto NEXT_MERGE_CHUNK;
    }
    finishworkingVars({});

    # finish up 
    printTimestampedMessage($message." done"); 
    close $inH_top;
    close $inH_bot;
    close $pileupH;
    close $snvH;
}

1;

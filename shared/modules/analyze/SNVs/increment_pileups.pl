use strict;
use warnings;

# functions that add aligments to stranded genome and read pileups
# threads alignments through a transient moving buffer for memory efficiency
# the following diagrams show how variant spans are mapped to the reference genome
#     matches and equal-length substitutions map one-to-one
#       --------------#------------  # indicates the reference position where variants are recorded
#       --------------:------------  :1
#       --------------*------------  12345*[actgn],[actgn]
#       --------------##-----------
#       --------------::-----------  :2 etc.
#       --------------**-----------  12345*[actgn][actgn],[actgn][actgn] etc.
#     simple deletions, 0 alt bases mapped to N ref bases
#       --------------##-----------
#       --------------DD-----------  12345*[actgn][actgn], etc.
#     indels with a net deletion, <N alt bases mapped to N ref bases
#       --------------##-----------  each # position reports variant *D
#       --------------*D-----------  12345*[actgn][actgn],[actgn] etc.
#     indels with a net insertion, >N aln bases mapped to N ref bases
#       --------------# -----------  each # position reports variant *I
#       --------------*I-----------  12345*[actgn],[actgn][actgn] etc.
#     simple insertions, N alt bases mapped to 0 ref bases, reported on the left-justified ref base
#       -------------#  -----------  # position reports adjacent variant II
#       --------------II-----------  12345*,[actgn][actgn] etc.
# based maps to read pileups are similar except that deletions, not insertions, are reported on the left-justified read base

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
    HIGH_QUAL_I     => 0,
    LOW_QUAL_I      => 1,
    HIGH_QUAL       => "Q",
    LOW_QUAL        => "q",
    #-------------
    CS_MATCH        => ":",
    CS_SUBSTITUTION => "*",
    CS_INSERTION    => "+",
    CS_DELETION     => "-",
    #-------------    
    NO_COVERAGE_RUN => "!q0Q0",  # used for no-coverage runs in initial stranded pileups (no-coverage runs omitted from final unstranded pileup)
    #-------------    
    CLIPPED         => "!",
    MASKED_LOW_QUAL => "q",
    #-------------
    READ1   => 1,
    READ2   => 2,
    #-------------
    TOP_STRAND      => 0, # fields shared between variants and variant coverage
    BOTTOM_STRAND   => 1, 
    #-------------
    N_HIGH_QUAL     => 2, # additional variant fields
    N_OBSERVED      => 3,
    REF_POS1        => 4,
    REF_BASES       => 5,
    ALT_BASES       => 6,
    COVERAGE        => 7,
    ZYGOSITY        => 8,
    GENOTYPE_COVERAGE_UNSTRANDED => 9,
    GENOTYPE_ZYGOSITY            => 10,
    #-------------
    UNSTRANDED      => 2, # additional coverage fields
    #-------------
    NA_ZYGOSITY     => 0, # zygosity codes
    SINGLETON       => 1,
    SUBCLONAL       => 2,
    HETEROZYGOUS    => 3,
    HOMOZYGOUS      => 4,
    NOT_IN_GENOTYPE => 5,
};

# variables
use vars qw(
    %outerEndpoints %workingVars %vars
    $pileupStartPos1  $lastProcessedI0 @pileup 
    $refCoverage_ $genotypeRefCoverage_ $MIN_READ_PILEUP_COVERAGE
    $chrom $chromSize $strandIndex0
    $pileupH $snvH
);
my $vaMatchOp_start = qr/^(\d+)(:)(\d+);/;
my $vaVarOp_start   = qr/^(\d+)(\*)([acgtn]*),([acgtn]*);/;
my $read1Map = [(READ1) x 1000]; # static maps for single-read alignments
my $read2Map = [(READ2) x 1000];
my ($readNMap, $invertRead2, $passedCoverage, $refPos1, $op, $opVal1, $opVal2);
my $nullAddedFields = join("\t", ("NA") x 3);
my $maskedBaseSpec = join("\t", MASKED_LOW_QUAL, MASKED_LOW_QUAL, $nullAddedFields);

# add the current grouped alignment to the stranded genome pileup buffer, starting from its leftmost strand position
sub incrementPileups_genome {
    my ($aln) = @_;
    my $leftPosI0   = $$aln[ALN_STRAND_POS1] - $pileupStartPos1;
    my $outerPos1_5 = $$aln[ALN_STRAND_POS1] - $$aln[ALN_READ_OFFET_5]; # identifying 5' most end of the read, including any clip
    my $vaTag = $$aln[ALN_VA_TAG];
    my @qsTag = split("", $$aln[ALN_QS_TAG]);
    while(
        $vaTag =~ s/$vaMatchOp_start// or
        $vaTag =~ s/$vaVarOp_start//
    ){
        my ($refPos1, $op, $opVal1, $opVal2) = ($1, $2, $3, $4);
        my $qualI = shift(@qsTag) ? HIGH_QUAL_I : LOW_QUAL_I; 
        if($op eq CS_MATCH){
            my $rightPosI0 = $leftPosI0 + $opVal1 - 1;
            foreach my $i0($leftPosI0..$rightPosI0){ # thus, record all values at all base positions in the match
                $pileup[$i0]{$op}[HIGH_QUAL_I] += $$aln[ALN_N_OBSERVED]; # refPos1 and nRefBases stripped, match is always :[HIGH_QUAL_I]
            }
            $leftPosI0 = $rightPosI0 + 1;
        } elsif(length($opVal1) > 0){ # all variants except isolated insertions
            my $varKey = $refPos1.$op.$opVal1.",".$opVal2;
            my $rightPosI0 = $leftPosI0 + length($opVal1) - 1;
            foreach my $i0($leftPosI0..$rightPosI0){ # again, same variant key and quality recorded at all applicable base positions in the variant stretch
                $pileup[$i0]{$varKey}[$qualI] += $$aln[ALN_N_OBSERVED]; # variant is more complex, e.g., 12345*ac,gta[HIGH_QUAL_I|LOW_QUAL_I]
            }
            $outerEndpoints{$varKey}{$outerPos1_5}++;
            $outerEndpoints{$outerPos1_5}{$varKey}++;
            $leftPosI0 = $rightPosI0 + 1;
        } else { # isolated insertions, no reference bases, always recorded on a single reference base
            # this is tricky since adding to a ref base that was already recorded
            # insertions recorded to left relative to the top strand (thus, right relative to the bottom strand index)
            my $varKey = $refPos1.$op.$opVal1.",".$opVal2;
            my $i0 = $leftPosI0 - ($strandIndex0 == TOP_STRAND ? 1 : 0);
            $pileup[$i0]{$varKey}[$qualI] += $$aln[ALN_N_OBSERVED];
            $outerEndpoints{$varKey}{$outerPos1_5}++;
            $outerEndpoints{$outerPos1_5}{$varKey}++;
        } 
    }
}

# add the current grouped alignment to the read-specific pileup buffer, starting from its 5'-most position
sub incrementPileups_read {
    my ($aln) = @_;

    # parse the base-level read number, which may change along the span of a merged read
    #   type          readN   xRead1   xRead2  from index_fragments.pl
    #   SE            1       1        0
    #   PE orphan     1       1        0
    #   PE merged     1       nR1>0    nR2>=0
    #   PE unmerged   1       1        0
    #   PE unmerged   2       0        1
    # ALN_X_READ1==TRUE means this is read1
    $readNMap = $$aln[ALN_X_READ1] ? 
        (   
            # ALN_X_READ1==1 means all aligned bases are from read1
            ($$aln[ALN_X_READ1] == 1 or $$aln[ALN_X_READ2] == 0) ?
                $read1Map :
                # o/w ALN_X_READ1 5' bases are read1, rest are read2
                [(READ1) x $$aln[ALN_X_READ1], (READ2) x $$aln[ALN_X_READ2]]
        ) : 
        # o/w this is read2 and all aligned bases are from read2
        $read2Map;

    # for merged reads, need to invert the read2 numbers
    $invertRead2 = $$readNMap[0] != $$readNMap[$#{$readNMap}];
    # 11111111122222
    # 1.......NN...1
    # TODO: some of this handling of read2 on merged reads may not be fully correct yet (rc, etc)

    # record any 5' end clips
    my $leftPosI0 = $$aln[ALN_READ_OFFET_5];
    if($leftPosI0){
        foreach my $i0(0..($leftPosI0 - 1)){ # thus, record all values at all base positions in the match
            addToReadPileup($aln, $i0, CLIPPED);
        }
    }

    # record aligned bases, include bases masked as low quality/unconfirmed variants
    my $vaTag = $$aln[ALN_VA_TAG];
    while(
        $vaTag =~ s/$vaMatchOp_start// or
        $vaTag =~ s/$vaVarOp_start//
    ){
        ($refPos1, $op, $opVal1, $opVal2) = ($1, $2, $3, $4);
        $passedCoverage = ( # restrict to bases with sufficient coverage
            counter_array_get($refCoverage_, $refPos1 - 1) >= $MIN_READ_PILEUP_COVERAGE and 
            (
                !$genotypeRefCoverage_ or # and, if requested, to bases with defined haplotypes in the genotype
                counter_array_get($genotypeRefCoverage_, $refPos1 - 1) == 2 # TODO: expose as EXPECTED_GENOTYPE_COVERAGE
            )
        );
        if($op eq CS_MATCH){
            my $rightPosI0 = $leftPosI0 + $opVal1 - 1;
            if($passedCoverage){
                foreach my $i0($leftPosI0..$rightPosI0){ # thus, record all values at all base positions in the match
                    addToReadPileup($aln, $i0, CS_MATCH);
                }
            }
            $leftPosI0 = $rightPosI0 + 1;
        } elsif(length($opVal2) > 0){ # all variants except isolated deletions
            my $rightPosI0 = $leftPosI0 + length($opVal2) - 1;
            if($passedCoverage){
                foreach my $i0($leftPosI0..$rightPosI0){ # again, same variant key and quality recorded at all applicable base positions in the variant stretch
                    addToReadPileup($aln, $i0);
                }
            }
            $leftPosI0 = $rightPosI0 + 1;
        } else { # isolated deletions, no read bases, always recorded on a single read base
            if($passedCoverage){
                my $i0 = $leftPosI0 - 1; # correct for read2 in merged read?
                addToReadPileup($aln, $i0);
            }
        } 
    }
}
sub addToReadPileup {
    my ($aln, $i0, $baseSpecKey) = @_;
    my $readN = $$readNMap[$i0];
    $invertRead2 and $readN == READ2 and $i0 = $#{$readNMap} - $i0;
    if($baseSpecKey){
        $baseSpecKey .= "\t".$baseSpecKey."\t".$nullAddedFields;
    } else {
        my $varKey = $refPos1.$op.$opVal1.",".$opVal2;
        my $var = $vars{$varKey};
        $baseSpecKey = $var ? 
            join("\t",
                $strandIndex0 == TOP_STRAND ? 
                    (      $opVal1,        $opVal2 ) : 
                    (getRc($opVal1), getRc($opVal2)), # correct for read2 in merged read?
                @{$var}[ZYGOSITY, GENOTYPE_ZYGOSITY, GENOTYPE_COVERAGE_UNSTRANDED]
            ) : 
            $maskedBaseSpec;
    }
    $pileup[$readN][$i0]{$baseSpecKey} += $$aln[ALN_N_OBSERVED];
}

# commit and purge a portion of the genome pileup buffer once no more reads are expected to overlap it
sub finishBases_genome {
    my ($lastPos1) = @_;
    my $lastI0 = $lastPos1 - $pileupStartPos1;
    my $hasBuffer = @pileup; # false on first process request on chromosome 
    my $workingRun = $hasBuffer ? $pileup[$lastProcessedI0] : ""; # the one run always present in the pileup buffer, that might be extended by the next processed position
    my $hasBreak;

    # process the collection of nominated baseSpec values into a final, streamlined pileup at each relevant reference base
    # every unique baseSpec must have the same set of quality counts at every base in the stretch, even when split into multiple runs, because:
    #   quality assessments for each read were made over all bases in the variant stretch (see getQsTag in parsing_functions.pl)
    #   all read observations of a variant stetch much have been encountered by the time any runs containing them are committed
    # however, we must still maintain independent low and high quality counts for each baseSpec, since we are yet to combine strands
    #   the other strand may have different qualities for the same variant stretch, potentially rescuing a low-quality variant on this strand
    foreach my $i0(($lastProcessedI0 + 1)..$lastI0){
        $pileup[$i0] = $pileup[$i0] ? join("|", map { 
            my $qN = $pileup[$i0]{$_};
            $_.LOW_QUAL.($$qN[LOW_QUAL_I] || 0).HIGH_QUAL.($$qN[HIGH_QUAL_I] || 0) # :q0Q55|12345*ac,gtaq11Q2 etc.
        } sort { $a cmp $b} keys %{$pileup[$i0]}) : NO_COVERAGE_RUN;        # or !q0Q0
        $pileup[$i0] ne $workingRun and $hasBreak = 1; # if never hits, we are just extending the current run
    }
    $lastProcessedI0 = $lastI0;
    ($hasBuffer and $hasBreak) or return;

    # collapse runs of identical baseSpec sets into single output lines
    TRY_NEXT_RUN: my ($runSpec, $runEndI0) = ($pileup[0], 0);
    while($runEndI0 < $lastProcessedI0 and $pileup[$runEndI0 + 1] eq $runSpec){
        $runEndI0++;
    }    
    if($runEndI0 < $lastProcessedI0){ # thus, always keep one run in the pileup buffer
        my $nPosInRun = $runEndI0 + 1;
        print $pileupH join("\t", $nPosInRun, $runSpec), "\n";
        splice(@pileup, 0, $nPosInRun);
        $pileupStartPos1 += $nPosInRun;
        $lastProcessedI0 -= $nPosInRun;
        goto TRY_NEXT_RUN;
    }
}

1;

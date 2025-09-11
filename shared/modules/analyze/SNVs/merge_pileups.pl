use strict;
use warnings;

# functions that merge strand-specific pileups into a single genome pileup

# constants
use constant {
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
    CS_MATCH        => ":",
    CS_SUBSTITUTION => "*",
    CS_INSERTION    => "+",
    CS_DELETION     => "-",
    #-------------
    NA_ZYGOSITY     => 0, # zygosity codes
    SINGLETON       => 1,
    SUBCLONAL       => 2,
    HETEROZYGOUS    => 3,
    HOMOZYGOUS      => 4,
    NOT_IN_GENOTYPE => 5,
    #-------------
    NA_COVERAGE     => 0, # coverage codes
    #-------------
    TOP_STRAND_COUNT    => "R",
    BOTTOM_STRAND_COUNT => "r",
};

# variables
use vars qw(
    $DEDUPLICATE_READS $IS_GENOTYPE_SNVS
    %outerEndpoints %workingVars %vars $genotypeRefCoverage_ $genotypeVarZygosity
    $pileupH $snvH $chromIndex1 $refCoverage_
);
my $baseSpec_all = qr/^(.+)q(\d+)Q(\d+)$/;
my $varKey_all = qr/(.+)\*(.*),(.*)/;
my $isolatedInsertion = qr/\*,/;
my $digit_start = qr/^\d/;

# pull one run of shared baseSpecs from either the top or bottom strand pileup
sub getStrandPileup {
    my ($inH, $refPos1) = @_;
    my $pileup = <$inH>;
    $pileup or return;
    chomp $pileup;
    my ($nPosInRun, $runSpec) = split("\t", $pileup);
    ($nPosInRun, $runSpec, $refPos1, $refPos1 + $nPosInRun - 1);
}

# add baseSpecs to the shared hash for the merging runs
sub fillBaseSpecs {
    my ($runSpec, $baseSpecs, $strandIndex0) = @_;
    foreach my $baseSpec(split('\|', $runSpec)){ # :q0Q55|12345*ac,gtaq11Q2 or !q0Q0
        my ($baseSpecKey, $nLowQual, $nHighQual) = ($baseSpec =~ m/$baseSpec_all/); # :q0Q55 or 12345*ac,gtaq11Q2 or !q0Q0
        $$baseSpecs{$baseSpecKey}[N_OBSERVED]    += $nLowQual + $nHighQual;
        $$baseSpecs{$baseSpecKey}[N_HIGH_QUAL]   += $nHighQual;
        $$baseSpecs{$baseSpecKey}[$strandIndex0] += $nLowQual + $nHighQual;
    }
}

# merge (portions of) runs from the top and bottom strand pileups
sub mergePileupRuns {
    my ($refPos1, $nPosInRun, $runSpec_top, $runSpec_bot) = @_;

    # merge the counts of all unique baseSpecs in the run
    my (%baseSpecs);
    fillBaseSpecs($runSpec_top, \%baseSpecs, TOP_STRAND);
    fillBaseSpecs($runSpec_bot, \%baseSpecs, BOTTOM_STRAND);

    # every unique baseSpec must have the same set of quality counts at every base in the stretch, even when split into multiple runs
    # thus, now that we have merged the strands, we can:
    #   require that all called variants have at least one high quality observation
    #   promote all observations of a variant to high quality if any observation of that variant is high quality
    #   reject all observations of a variant if no observation of that variant is high quality
    #   rejected baseSpecs, i.e., varKeys, are:
    #      ignored as artifacts
    #      omitted from the final pileup as if those bases never existed
    my @coverage = (0, 0, 0);
    foreach my $baseSpecKey(keys %baseSpecs){

        # record the coverage of accepted baseSpecs for VAF calculation
        # isolated insertions are:
        #     identifed by string *, in a baseSpecKey i.e., no reference bases
        #     not counted in coverage as they were recorded on the same reference base as an adjacent match
        # all other variants including deletions and indels increase coverage 
        #     of all reference bases in the variant stretch (no match was recorded here)
        #     whether or not the variant is the same length as the reference span
        if($baseSpecs{$baseSpecKey}[N_HIGH_QUAL]){ # TODO: expose the required number/frequency of high quality observations as an option
            unless($baseSpecKey =~ m/$isolatedInsertion/){
                $coverage[TOP_STRAND]    += ($baseSpecs{$baseSpecKey}[TOP_STRAND]    || 0); 
                $coverage[BOTTOM_STRAND] += ($baseSpecs{$baseSpecKey}[BOTTOM_STRAND] || 0); 
                $coverage[UNSTRANDED]    +=  $baseSpecs{$baseSpecKey}[N_OBSERVED];
            }

            # record the leftmost position of accepted variant stretches
            $baseSpecKey =~ m/$digit_start/ and 
                @{$baseSpecs{$baseSpecKey}}[REF_POS1, REF_BASES, ALT_BASES] = ($baseSpecKey =~ m/$varKey_all/);
            
        # permanently remove rejected baseSpecs with no high quality observations
        } else {
            delete $baseSpecs{$baseSpecKey};  
            $outerEndpoints{$baseSpecKey} and delete $outerEndpoints{$baseSpecKey};
        }
    }
    $coverage[UNSTRANDED] or return $refPos1 + $nPosInRun; # uncovered stretches are not included in the final pileup (unlike the interim pileup above)
    counter_array_set_range($refCoverage_, $refPos1 - 1, $nPosInRun, $coverage[UNSTRANDED]);

    # commit the unstranded pileup, merged over both strands
    # accumulate working SNVs prior to committing them once we are past their span
    # let the run with the highest coverage describe SNVs that are split over multiple runs
    print $pileupH join("\t", 
        $chromIndex1, $refPos1 - 1, $refPos1 - 1 + $nPosInRun, # BED format coordinates of the run
        join("|", map {
            my $bs = $baseSpecs{$_};

            # accumate SNV metadata
            if($_ ne CS_MATCH){ # only 12345*ac,gta remains
                if($workingVars{$_}){
                    if($coverage[UNSTRANDED] > $workingVars{$_}[COVERAGE][UNSTRANDED]){
                        $$bs[COVERAGE] = \@coverage;
                        $workingVars{$_} = $bs;
                    }
                } else {
                    $$bs[COVERAGE] = \@coverage;
                    $workingVars{$_} = $bs; # thus has keys coverage, nObserved, nHighQual
                }
            }

            # return the merged runSpec being printed
            # retains strand counts only, all variants have been validated by high quality observations
            # :R35r20|12345*ac,gtaR11r4
            # due to assembly errors, alignment errors, or copy number differences, `genotype SNVs` coverage can exceed two
            # here is an example from minimap2 alignment of an NA12878 assembly to hs1
            # three aligned contigs cover a region from hs1 chr1 base 233905 to 246299
            #   -       chr1    218835  246299  27320   27495   60 # from minimap2 PAF
            #   +       chr1    233019  859628  621115  631407  60
            #   +       chr1    233904  859248  619408  628345  60
            #   1       233857  233904  :R1r1 # from genome pileup
            #   1       233904  233943  :R2r1
            #   ...
            #   1       246220  246299  :R2r1
            #   1       246299  246347  :R2r0
            $_.TOP_STRAND_COUNT.($$bs[TOP_STRAND] || 0).BOTTOM_STRAND_COUNT.($$bs[BOTTOM_STRAND] || 0);
        } sort { $a cmp $b} keys %baseSpecs)
    ), "\n";

    # return the next position to process
    finishworkingVars(\%baseSpecs);
    $refPos1 + $nPosInRun;
}

# finish SNVs that are no longer being observed, i.e., that we have moved past
sub finishworkingVars {
    my ($pendingVars) = @_;
    foreach my $varKey(keys %workingVars){
        $$pendingVars{$varKey} and next; # this SNV stretch has not necessarily ended yet, still working, wait

        # use the outerEndpoints hash to potentially drop SNVs that appear to arise by PCR amplification
        # at present we don't attempt to reduce counts of reported SNVS, only to drop inhomogenous SNVs that are not validated by a second outer endpoint
        #   $outerEndpoints{$varKey}{$outerPos1_5}++; # the collection of outer endpoints for each SNV stretch
        #   $outerEndpoints{$outerPos1_5}{$varKey}++; # the collection of SNV stretches for each outer endpoint
        # --------------------------
        # -----------X--------------
        #      ------X-------------------
        #      ------X-------------------
        if($DEDUPLICATE_READS and keys %{$outerEndpoints{$varKey}} == 1){ # if a second outer endpoint validates the SNV, keep all counts, nothing to do
            my $outerPos1_5 = (keys %{$outerEndpoints{$varKey}})[0];
            if($outerEndpoints{$outerPos1_5}{$varKey} / $workingVars{$varKey}[COVERAGE][UNSTRANDED] <= 0.5){ # all or most reads reported the SNV, keep it as is
                delete $workingVars{$varKey}; # a single outer endpoint reported the SNV, but a minority of those reads reported the SNV, drop it as untrustworthy
                delete $outerEndpoints{$varKey};
                next;
            }
        }

        # finish describing the SNV stretch
        if($IS_GENOTYPE_SNVS){
            $workingVars{$varKey}[ZYGOSITY] = $workingVars{$varKey}[N_OBSERVED] == $workingVars{$varKey}[COVERAGE][UNSTRANDED] ? 
                HOMOZYGOUS :
                HETEROZYGOUS;
        } else {
            my $vaf = $workingVars{$varKey}[N_OBSERVED] / $workingVars{$varKey}[COVERAGE][UNSTRANDED];
            $workingVars{$varKey}[ZYGOSITY] = $workingVars{$varKey}[N_OBSERVED] == 1 ? # this is crude, not intended to be perfect
                SINGLETON :                                                            # in particular, VAF here does not account for coverage and therefore confidence
                ( 
                    $vaf > 0.75 ? 
                        HOMOZYGOUS: 
                        (
                            $vaf > 0.25 ? 
                                HETEROZYGOUS : 
                                SUBCLONAL
                        )
                );
                $workingVars{$varKey}[GENOTYPE_COVERAGE_UNSTRANDED] = $genotypeRefCoverage_ ?
                    counter_array_get($genotypeRefCoverage_, $workingVars{$varKey}[REF_POS1] - 1) :
                    NA_COVERAGE;
                $workingVars{$varKey}[GENOTYPE_ZYGOSITY] = $genotypeVarZygosity ?
                    ($$genotypeVarZygosity{$varKey} || NOT_IN_GENOTYPE) :
                    NA_ZYGOSITY;
        }

        # commit the SNV stretch and clean up
        $vars{$varKey} = $workingVars{$varKey};
        delete $workingVars{$varKey}; # done with this SNV stretch, it cannot be encountered again
        delete $outerEndpoints{$varKey};
        # TODO: change this to VCF?
        print $snvH join("\t", 
            $chromIndex1,
            $vars{$varKey}[REF_POS1] - 1, 
            $vars{$varKey}[REF_POS1] - 1 + max(1, length($vars{$varKey}[REF_BASES])), # BED format coordinates of the SNV
            $vars{$varKey}[REF_BASES], 
            $vars{$varKey}[ALT_BASES], 
            (map { $vars{$varKey}[$_] || 0 } TOP_STRAND, BOTTOM_STRAND, N_HIGH_QUAL, N_OBSERVED),
            (map { $vars{$varKey}[COVERAGE][$_] } TOP_STRAND, BOTTOM_STRAND, UNSTRANDED),
            $vars{$varKey}[ZYGOSITY],
            ($IS_GENOTYPE_SNVS ? () : map { $vars{$varKey}[$_] } GENOTYPE_COVERAGE_UNSTRANDED, GENOTYPE_ZYGOSITY)
        ), "\n";
    }
}

1;

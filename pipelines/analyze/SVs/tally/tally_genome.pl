use strict;
use warnings;

# actions:
#   aggregate genome spans split by target region and gene
# input:
#   piped from tally_genome.sh
# output:
#   various counts in a log report

# initialize reporting
our $script = "tally_genome";
our $error  = "$script error";
my (%nBases, %totalBases);

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
resetCountFile();

# constants
use constant {
    BED_CHROM      => 0,
    BED_START0     => 1,
    BED_END1       => 2,
    BED_NAME       => 3,
    BED_SCORE      => 4,
    BED_STRAND     => 5,
    #-------------
    TRUE_STR  => 'true',
    FALSE_STR => 'false',
    NA_STR    => 'NA',
};

# run spans
while (my $span = <STDIN>){
    chomp $span;
    my @span = split("\t", $span);

    # set the single genome, or one of two mixed genomes
    my ($chrom, $genome) = split("_", $span[BED_CHROM]);
    $genome or $genome = $ENV{GENOME};

    # parse the match to genome segments statified by target region and gene
    # "gap" refers to regions that are not in a gene or a target region
    my $nBases = $span[BED_END1] - $span[BED_START0];
    my %parts    = (TARGET => {}, GENE => {}, GAP => {});
    my %hasParts = (TARGET => 0,  GENE => 0,  GAP => 0);
    map {
        my ($type, $name) = split("_", $_);
        $parts{$type}{$name || "GAP"}++;
        $hasParts{$type}++;
    } split(",", $span[BED_NAME]);
    my $isGap    = $hasParts{GAP}    ? TRUE_STR : FALSE_STR;
    my $onTarget = $hasParts{TARGET} ? TRUE_STR : FALSE_STR;
    my $inGene   = $hasParts{GENE}   ? TRUE_STR : FALSE_STR;
    my $targetGene = $hasParts{TARGET} ? (keys %{$parts{TARGET}})[0] : NA_STR;
    my $inTargetGene = ($hasParts{TARGET} and $hasParts{GENE}) ? (
        $parts{GENE}{$targetGene} ? TRUE_STR : FALSE_STR
    ) : FALSE_STR;

    # count the bases in each type of region
    my $key = join("\t", 
        $genome, $isGap, $onTarget, $inGene, $inTargetGene, $targetGene
    );
    $nBases{$key} += $nBases;
    $totalBases{$genome} += $nBases;
    $totalBases{ALL} += $nBases;
}

# print tally
print join("\t", qw(
    genome isGap onTarget inGene inTargetGene targetGene 
    nBases fracBases_genome fracBases_all
)), "\n";
foreach my $key (sort { $b cmp $a } keys %nBases){
    my ($genome) = split("\t", $key);
    my $fracBases_genome = int($nBases{$key} / $totalBases{$genome} * 10000) / 10000;
    my $fracBases_all    = int($nBases{$key} / $totalBases{ALL}     * 10000) / 10000;
    print join("\t",
        $key, 
        $nBases{$key}, $fracBases_genome, $fracBases_all
    ), "\n";
}

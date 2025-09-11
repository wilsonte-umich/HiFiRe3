use strict;
use warnings;

# actions:
#   report some useful summary information to answer questions like:
#       what is the rate of SV junction artifacts per base?
# input:
#   $SV_ALIGNMENTS_FILE
#   $SV_FINAL_JUNCTIONS_FILE_1
# output:
#   various counts in a log report

# initialize reporting
our $script = "tally_junctions";
our $error  = "$script error";
my ($nUniqJxns, $nJxnsN3Plus, $nSingletonTrans, $nInterGenome,
    $nOnTarget1, $nOnTarget2, $nInGene1, $nInGene2,
    $nInTargetGene1, $nInTargetGene2) = (0) x 100;
my (%jxnTypes, %nReads, %nBases);
my (%nBkpts, %totalBkpts);

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
resetCountFile();

# environment variables
fillEnvVar(\our $SV_ALIGNMENTS_FILE,        'SV_ALIGNMENTS_FILE');
fillEnvVar(\our $SV_FINAL_JUNCTIONS_FILE_1, 'SV_FINAL_JUNCTIONS_FILE_1');
fillEnvVar(\our $GENOME,                    'GENOME');
fillEnvVar(\our $IS_COMPOSITE_GENOME,       'IS_COMPOSITE_GENOME');

# constants
use constant {
    ALN_CHROM_INDEX1_1      => 0,
    ALN_REF_POS1_5          => 1,
    ALN_REF_POS1_2          => 2,
    ALN_REF_POS1_3          => 3,
    ALN_STRAND_INDEX1_0     => 4,
    ALN_JXN_TYPE            => 5,
    ALN_HAPLOTYPES          => 6,
    ALN_PATH_N              => 7,
    ALN_N_JUNCTIONS         => 8,
    ALN_ALN_N               => 9,
    ALN_PATH_READ_HAS_SV    => 10,
    ALN_N_OBSERVED          => 11,
    #-------------
    OJXN_CHROM_INDEX1_1     => 0, # input columns
    OJXN_REF_POS1_1         => 1,
    OJXN_STRAND_INDEX0_1    => 2,
    OJXN_CHROM_INDEX1_2     => 3,
    OJXN_REF_POS1_2         => 4,
    OJXN_STRAND_INDEX0_2    => 5,
    OJXN_JXN_TYPE           => 6,
    UJXN_N_OBSERVED         => 7,
    UJXN_ALN_OFFSET         => 8,
    UJXN_JXN_BASES          => 9,
    UJXN_PATHS              => 10,
    UJXN_N_PATH_JUNCTIONS   => 11,
    UJXN_READ_HAS_SV        => 12,
    UJXN_JSRC_MAPQ          => 13, 
    UJXN_JSRC_DE_TAG        => 14,
    UJXN_JSRC_SITE_DIST     => 15,
    UJXN_JSRC_QNAMES        => 16, 
    UJXN_JSRC_SEQS          => 17,
    UJXN_JSRC_QUALS         => 18,
    UJXN_JSRC_CIGARS        => 19,
    UJXN_JSRC_ORIENTATIONS  => 20,
    UJXN_HAS_ALT_ALIGNMENT  => 21,
    UJXN_SV_SIZE            => 22,
    UJXN_IS_INTERGENOME     => 23,
    UJXN_TARGET_1           => 24,
    UJXN_TARGET_DIST_1      => 25,
    UJXN_TARGET_2           => 26,
    UJXN_TARGET_DIST_2      => 27,
    UJXN_GENES_1            => 28,
    UJXN_GENE_DIST_1        => 29,
    UJXN_GENES_2            => 30,
    UJXN_GENE_DIST_2        => 31,
    UJXN_IS_EXCLUDED_1      => 32,
    UJXN_IS_EXCLUDED_2      => 33,
    UJXN_BKPT_COVERAGE_1    => 34,
    UJXN_BKPT_COVERAGE_2    => 35,
    #-------------
    EXPECTED   => 0, # junction sequencing orientations
    UNEXPECTED => 1,
    #-------------
    TYPE_PROPER         => 0, # junction/edge types
    TYPE_DELETION       => 1,
    TYPE_DUPLICATION    => 2,
    TYPE_INVERSION      => 3,
    TYPE_TRANSLOCATION  => 4,
    #-------------
    TRUE_STR  => 'true',
    FALSE_STR => 'false',
    NULL_STR  => '*',
    NA_STR    => 'NA',
};
my @jxnTypeNames = qw(proper deletion duplication inversion translocation);

# initialize the genome
use vars qw(%revChromIndex);
setCanonicalChroms();

# tally junctions for total counts in log file
open my $jxnH, "-|", "zcat $SV_FINAL_JUNCTIONS_FILE_1" or die "could not open $SV_FINAL_JUNCTIONS_FILE_1: $!\n";
while (my $jxn = <$jxnH>){
    chomp $jxn;
    my @jxn = split("\t", $jxn);

    # junction counts by validation status
    $nUniqJxns++;
    $jxn[UJXN_N_OBSERVED] >= 3 and $nJxnsN3Plus++;
    $jxn[UJXN_N_OBSERVED] == 1 and $jxn[OJXN_JXN_TYPE] == TYPE_TRANSLOCATION and $nSingletonTrans++;

    # count of inter-genome unequivocal artifacts
    $jxn[UJXN_IS_INTERGENOME] and $nInterGenome++;

    # junction counts by target/gene match
    ($jxn[UJXN_TARGET_1] ne NULL_STR or  $jxn[UJXN_TARGET_2] ne NULL_STR) and $nOnTarget1++;
    ($jxn[UJXN_TARGET_1] ne NULL_STR and $jxn[UJXN_TARGET_2] ne NULL_STR) and $nOnTarget2++;
    ($jxn[UJXN_GENES_1] ne ","  or  $jxn[UJXN_GENES_2] ne ",")  and $nInGene1++;
    ($jxn[UJXN_GENES_1] ne ","  and $jxn[UJXN_GENES_2] ne ",")  and $nInGene2++;
    ($jxn[UJXN_TARGET_1] ne NULL_STR or  $jxn[UJXN_TARGET_2] ne NULL_STR) and
    ($jxn[UJXN_GENES_1] ne ","  or  $jxn[UJXN_GENES_2] ne ",")  and 
    $nInTargetGene1++;
    ($jxn[UJXN_TARGET_1] ne NULL_STR and $jxn[UJXN_TARGET_2] ne NULL_STR) and
    ($jxn[UJXN_GENES_1] ne ","  and $jxn[UJXN_GENES_2] ne ",")  and 
    $nInTargetGene2++;

    # tally junctions for more expanded analysis in tally file
    # importantly, not all junctions are tallied in this section by design!
    if(
        ($jxn[UJXN_N_OBSERVED] == 1 or $jxn[UJXN_N_OBSERVED] >= 3) and # do not tally jxns with exactly 2 reads
        ($jxn[OJXN_JXN_TYPE] == TYPE_TRANSLOCATION or $jxn[UJXN_SV_SIZE] > 1000) # ignore super-small SVs
    ){

        # determine the genome of each breakpoint
        my ($chrom1, $genome1) = split("_", $revChromIndex{$jxn[OJXN_CHROM_INDEX1_1]});
        $genome1 or $genome1 = $GENOME;
        my ($chrom2, $genome2) = split("_", $revChromIndex{$jxn[OJXN_CHROM_INDEX1_2]});
        $genome2 or $genome2 = $GENOME;

        # parse the match to genome segments stratified by target region and gene
        my $onTarget1 = $jxn[UJXN_TARGET_1] ne NULL_STR ? TRUE_STR : FALSE_STR;
        my $inGene1   = $jxn[UJXN_GENES_1]  ne "," ? TRUE_STR : FALSE_STR;
        my $isGap1    = ($jxn[UJXN_TARGET_1] eq NULL_STR and $jxn[UJXN_GENES_1] eq ",") ? TRUE_STR : FALSE_STR;
        my $onTarget2 = $jxn[UJXN_TARGET_2] ne NULL_STR ? TRUE_STR : FALSE_STR;
        my $inGene2   = $jxn[UJXN_GENES_2]  ne "," ? TRUE_STR : FALSE_STR;
        my $isGap2    = ($jxn[UJXN_TARGET_2] eq NULL_STR and $jxn[UJXN_GENES_2] eq ",") ? TRUE_STR : FALSE_STR;
        my $uxnTarget1 = $jxn[UJXN_TARGET_1]; # must assign new variables to allow proper regex parsing
        my $uxnTarget2 = $jxn[UJXN_TARGET_2];
        my $inTargetGene1 = ($uxnTarget1 ne NULL_STR and $jxn[UJXN_GENES_1] ne ",") ? (
            $jxn[UJXN_GENES_1] =~ m/,$uxnTarget1,/ ? TRUE_STR : FALSE_STR
        ) : FALSE_STR;
        my $inTargetGene2 = ($uxnTarget2 ne NULL_STR and $jxn[UJXN_GENES_2] ne ",") ? (
            $jxn[UJXN_GENES_2] =~ m/,$uxnTarget2,/ ? TRUE_STR : FALSE_STR
        ) : FALSE_STR;

        # count breakpoints by junction type and validation status
        my $jxnType       = $jxnTypeNames[$jxn[OJXN_JXN_TYPE]];
        my $isLocalSv     = ($jxn[OJXN_JXN_TYPE] != TYPE_TRANSLOCATION and $jxn[UJXN_SV_SIZE] < 1e6) ? TRUE_STR : FALSE_STR;
        my $isInterGenome = $jxn[UJXN_IS_INTERGENOME] ? TRUE_STR : FALSE_STR;
        my $isValidated   = $jxn[UJXN_N_OBSERVED] == 1 ? FALSE_STR : TRUE_STR; # if not 1, then 3 or more per above
        my $isExpected    = $jxn[UJXN_READ_HAS_SV] ? FALSE_STR : TRUE_STR; # may be false inappropriately for multi-junction reads

        # assemble breakpoint category keys and count
        my $key1 = join("\t", 
            $genome1, $isGap1, $onTarget1, $inGene1, $inTargetGene1, 
            $jxn[UJXN_TARGET_1] eq NULL_STR ? NA_STR : $jxn[UJXN_TARGET_1],
            $jxnType, $isLocalSv, $isInterGenome, $isValidated, $isExpected
        );
        my $key2 = join("\t", 
            $genome2, $isGap2, $onTarget2, $inGene2, $inTargetGene2, 
            $jxn[UJXN_TARGET_2] eq NULL_STR ? NA_STR : $jxn[UJXN_TARGET_2],
            $jxnType, $isLocalSv, $isInterGenome, $isValidated, $isExpected
        );
        $nBkpts{$key1}++;
        $nBkpts{$key2}++;
        $totalBkpts{$genome1}++;
        $totalBkpts{$genome2}++;
        $totalBkpts{ALL}++;
    }
}
close $jxnH;

# print tally to STDOUT
print join("\t", qw(
    genome isGap onTarget inGene inTargetGene 
    targetGene 
    jxnType isLocalSv isInterGenome isValidated isExpected
    nBkpts fracBkpts_genome fracBkpts_all
)), "\n";
foreach my $key (sort { $b cmp $a } keys %nBkpts){
    my ($genome) = split("\t", $key);
    my $fracBkpts_genome = int($nBkpts{$key} / $totalBkpts{$genome} * 10000) / 10000;
    my $fracBkpts_all    = int($nBkpts{$key} / $totalBkpts{ALL}     * 10000) / 10000;
    print join("\t",
        $key, 
        $nBkpts{$key}, $fracBkpts_genome, $fracBkpts_all
    ), "\n";
}

# print summary information
printCount(commify($nUniqJxns),         'nUniqJxns',        'final unique junctions');
printCount(commify($nJxnsN3Plus),       'nJxnsN3Plus',      'confirmed unique junctions, nObserved >= 3');
printCount(commify($nSingletonTrans),   'nSingletonTrans',  'singleton translocations, presumed artifacts');
$IS_COMPOSITE_GENOME and 
printCount(commify($nInterGenome),      'nInterGenome',     'inter-genome translocations, unambiguous artifacts');

printCount(commify($nOnTarget1),        'nOnTarget1',       'junctions on target, at least one breakpoint');
printCount(commify($nOnTarget2),        'nOnTarget2',       'junctions on target, both breakpoints');
printCount(commify($nInGene1),          'nInGene1',         'junctions in gene, at least one breakpoint');
printCount(commify($nInGene2),          'nInGene2',         'junctions in gene, both breakpoints');

printCount(commify($nInTargetGene1),    'nInTargetGene1',   'junctions in target gene, at least one breakpoint');
printCount(commify($nInTargetGene2),    'nInTargetGene2',   'junctions in target gene, both breakpoints');

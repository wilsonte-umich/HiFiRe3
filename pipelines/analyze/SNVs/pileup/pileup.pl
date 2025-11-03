use strict;
use warnings;

# actions:
#   combine indexed single alignments into genome pileup
#   assess the presence of (sub)clonal (and thus expected) SNVs/indel alleles
#   if requested, compare observed SNV
#   work is parallelized by chromosome
# input:
#   $SNV_ALNS_PREFIX.<zeroPaddedChromIndex1>.txt.gz
# output:
#   files per chromosome, for subsequent aggregation:
#       SNV_GENOME_PILEUP_PREFIX.* = unstranded merged pileup bed.bgz file
#       SNV_SUMMARY_TABLE_PREFIX.* = called SNV/indels bed.bgz file, including singletons and (sub)clonal alleles, optionally with expected genotypes

# initialize reporting
our $script = "pileup";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric counter_array);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);

# environment variables
fillEnvVar(\our $N_CPU,                     'N_CPU');
fillEnvVar(\our $GENOME_FASTA,              'GENOME_FASTA');
fillEnvVar(\our $SNV_ALNS_PREFIX,           'SNV_ALNS_PREFIX');
fillEnvVar(\our $TMP_PILEUP_DIR,            'TMP_PILEUP_DIR');
fillEnvVar(\our $SNV_GENOME_PILEUP_PREFIX,  'SNV_GENOME_PILEUP_PREFIX');
fillEnvVar(\our $SNV_SUMMARY_TABLE_PREFIX,  'SNV_SUMMARY_TABLE_PREFIX');
fillEnvVar(\our $ACTION_DIR,                'ACTION_DIR');
fillEnvVar(\our $MODULES_DIR,               'MODULES_DIR');
fillEnvVar(\our $MIN_READ_PILEUP_COVERAGE,  'MIN_READ_PILEUP_COVERAGE');
fillEnvVar(\our $GENOTYPE_SNV_BASE_COVERAGE_PREFIX, 'GENOTYPE_SNV_BASE_COVERAGE_PREFIX', 1, "");
fillEnvVar(\our $GENOTYPE_SNV_SUMMARY_TABLE,        'GENOTYPE_SNV_SUMMARY_TABLE', 1, "");

# constants
use constant {
    TOP_STRAND      => 0, # fields shared between variants and variant coverage
    BOTTOM_STRAND   => 1, 
    #-------------    
    READ1   => 1,
    READ2   => 2,
    #-------------
    SNV_CHROM_INDEX1                 => 0,  # base positions on chrom affected by SNV/indel, BED format
    SNV_START0                       => 1,
    SNV_END1                         => 2,
    SNV_REF_BASES                    => 3,  # ref and alt base strings defining the SNV/indel
    SNV_ALT_BASES                    => 4,
    SNV_N_TOP_STRAND                 => 5,  # observation counts, i.e., number of reads supporting the SNV/indel
    SNV_N_BOTTOM_STRAND              => 6,
    SNV_N_HIGH_QUAL                  => 7,
    SNV_N_OBSERVED                   => 8,
    SNV_COVERAGE_TOP_STRAND          => 9,  # number of reads overlapping the SNV/indel, whether or not they support it
    SNV_COVERAGE_BOTTOM_STRAND       => 10,
    SNV_COVERAGE_UNSTRANDED          => 11,
    SNV_ZYGOSITY                     => 12, # inferred genotype/mosacism derived from SNV_N_OBSERVED and SNV_COVERAGE_UNSTRANDED
    SNV_GENOTYPE_COVERAGE_UNSTRANDED => 13, # coverage and zygosity from external genotype information (absent from genotype SNV tables)
    SNV_GENOTYPE_ZYGOSITY            => 14,
};

# initialize the genome
print STDERR "collecting chrom data\n";
use vars qw(%chromIndex);
setCanonicalChroms();
our @nuclearChroms = getNuclearChroms();
our %chromSizes = getChromSizes("$GENOME_FASTA.fai");

# parallelize by chromosome for speed
# user responsible for setting N_CPU to prevent memory overruns
print STDERR "piling up read alignments by chromosome\n";
$| = 1;
launchChildThreads(\&processChroms);
use vars qw(@readH @writeH);
my $writeH = $writeH[1];
# foreach my $chrom("chr21", "chr22"){
foreach my $chrom(@nuclearChroms){
    my $chromIndex1 = $chromIndex{$chrom} or next;
    my $chromSize = $chromSizes{$chrom} or die "$error: could not find size for $chrom\n";
    $writeH = $writeH[$chromIndex1 % $N_CPU + 1];
    print $writeH join("\t", $chrom, $chromIndex1, $chromSize), "\n";
}
finishChildThreads();

# child process to parse bam read pairs
our ($chrom, $chromIndex1, $paddedChromIndex1, $chromSize, $strandIndex0, $childN_); # reset by thread line
our (%outerEndpoints, %workingVars, %vars, $genotypeVarZygosity); # reset every chrom
our ($pileupStartPos1, $lastProcessedI0, @pileup, $refCoverage_, $genotypeRefCoverage_); # reset as needed for pileup context
our ($pileupH, $snvH);
sub processChroms {
    my ($childN) = @_;
    my $readH = $readH[$childN];
    map { require "$MODULES_DIR/analyze/SNVs/$_.pl" } qw(genome_pileups increment_pileups merge_pileups);

    # process data by chrom
    while(my $line = <$readH>){
        chomp $line;
        ($chrom, $chromIndex1, $chromSize) = split("\t", $line);
        $paddedChromIndex1 = sprintf("%02d", $chromIndex1);
        (%outerEndpoints, %workingVars, %vars) = (); # collected across both strands per chrom

        # pileup is only on one effective strand for duplex hifi reads

        runStrandedGenomePileup($childN, TOP_STRAND,    \&getSnvAlnsStream);
        runStrandedGenomePileup($childN, BOTTOM_STRAND, \&getSnvAlnsStream);

        $refCoverage_ = counter_array_init($chromSize);
        $genotypeRefCoverage_ = $GENOTYPE_SNV_BASE_COVERAGE_PREFIX ? counter_array_init_from_file(
            "$GENOTYPE_SNV_BASE_COVERAGE_PREFIX.$paddedChromIndex1.raw"
        ) : undef;
        $genotypeVarZygosity = loadGenotypeVarZygosity();
        mergeStrandedGenomePileups($childN);


        # $genotypeVarZygosity = undef;
        # @pileup = ();
        # runReadPileup($childN, TOP_STRAND);
        # runReadPileup($childN, BOTTOM_STRAND);
        # counter_array_release($refCoverage_, 1);
        # $genotypeRefCoverage_ and counter_array_release($genotypeRefCoverage_, 1);
        # printReadPileup();
    }
}
sub getSnvAlnsStream {
    "zcat $SNV_ALNS_PREFIX.$paddedChromIndex1.$strandIndex0.txt.gz";
}
sub loadGenotypeVarZygosity {
    my $genotypeVarZygosity = {};
    $GENOTYPE_SNV_SUMMARY_TABLE or return $genotypeVarZygosity;
    open my $genotypeH, "-|", "tabix $GENOTYPE_SNV_SUMMARY_TABLE $chromIndex1"
        or die "$error: could not open: $GENOTYPE_SNV_SUMMARY_TABLE: $!\n";
    while(my $line = <$genotypeH>){
        chomp $line;
        my @gtVar = split("\t", $line);
        my $varKey = ($gtVar[SNV_START0] + 1)."*$gtVar[SNV_REF_BASES],$gtVar[SNV_ALT_BASES]";
        $$genotypeVarZygosity{$varKey} = $gtVar[SNV_ZYGOSITY];
    }
    close $genotypeH;
    $genotypeVarZygosity;
}

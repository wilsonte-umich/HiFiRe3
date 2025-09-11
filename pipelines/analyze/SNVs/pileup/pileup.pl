use strict;
use warnings;

# actions:
#   combine indexed 5'-most alignments into two genome pileups, one for each strand (temporary files only used by this script)
#   merged the two strand pileups into a single genome pileup (only this merged genome pileup is retained as script output)
#   use the above information to combine indexed 5'-most alignments into one or two additional read pileups
#   assess the presence of (sub)clonal (and thus expected) SNVs/indel alleles
#   if requested, compare observed SNV
#   work is parallelized by chromosome
# input:
#   $SNV_ALNS_PREFIX.<zeroPaddedChromIndex1>.<strandIndex0>.txt.gz
# output:
#   three files per chromosome, for subsequent aggregation:
#       SNV_GENOME_PILEUP_PREFIX.* = unstranded merged pileup bed.bgz file
#       SNV_SUMMARY_TABLE_PREFIX.* = called SNV/indels bed.bgz file, including singletons and (sub)clonal alleles, optionally with expected genotypes
#       SNV_READ_PILEUP_PREFIX.*  = read-specific pileup txt.gz file (read1 plus read2 when applicable for paired-end data)

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
fillEnvVar(\our $SNV_READ_PILEUP_PREFIX,    'SNV_READ_PILEUP_PREFIX');
fillEnvVar(\our $SNV_SUMMARY_TABLE_PREFIX,  'SNV_SUMMARY_TABLE_PREFIX');
fillEnvVar(\our $DEDUPLICATE_READS,         'DEDUPLICATE_READS');
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
        runStrandedGenomePileup($childN, TOP_STRAND,    \&getSnvAlnsStream);
        runStrandedGenomePileup($childN, BOTTOM_STRAND, \&getSnvAlnsStream);
        $refCoverage_ = counter_array_init($chromSize);
        $genotypeRefCoverage_ = $GENOTYPE_SNV_BASE_COVERAGE_PREFIX ? counter_array_init_from_file(
            "$GENOTYPE_SNV_BASE_COVERAGE_PREFIX.$paddedChromIndex1.raw"
        ) : undef;
        $genotypeVarZygosity = loadGenotypeVarZygosity();
        mergeStrandedGenomePileups($childN);
        $genotypeVarZygosity = undef;
        @pileup = ();
        runReadPileup($childN, TOP_STRAND);
        runReadPileup($childN, BOTTOM_STRAND);
        counter_array_release($refCoverage_, 1);
        $genotypeRefCoverage_ and counter_array_release($genotypeRefCoverage_, 1);
        printReadPileup();
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

# run the indexed alignment a second time to create a read pileup
# use the catalog of (sub)clonal SNVs to assign variants as expected or unexpected
sub runReadPileup {
    ($childN_, $strandIndex0) = @_;
    my $message = "thread $childN_\trunReadPileup $chrom $strandIndex0";
    printTimestampedMessage($message); 

    # open file handles
    my $snvAlnsFile = "$SNV_ALNS_PREFIX.$paddedChromIndex1.$strandIndex0.txt.gz";
    open my $snvAlnsH, "-|", "zcat $snvAlnsFile" 
        or die "$error: could not open: $snvAlnsFile: $!\n";

    # run the alignments
    while(my $line = <$snvAlnsH>){
        chomp $line;
        my @aln = split("\t", $line);
        incrementPileups_read(\@aln);
    }

    # finish up 
    printTimestampedMessage($message." done"); 
    close $snvAlnsH;
}

# print the combined strand read pileups
sub printReadPileup {
    open $pileupH, "|-", "gzip -c > $SNV_READ_PILEUP_PREFIX.$paddedChromIndex1.txt.gz"
        or die "$error: could not open: $!\n";
    foreach my $readN(READ1, READ2){
        $pileup[$readN] or next;
        foreach my $i0(0..$#{$pileup[$readN]}){
            my $pileup = $pileup[$readN][$i0];
            $pileup or next;
            foreach my $baseSpecKey(keys %$pileup){
                print $pileupH join("\t", $readN, $i0 + 1, $baseSpecKey, $$pileup{$baseSpecKey}), "\n";
            }
        }
    }
    close $pileupH;
}


# [wilsonte@gl-login2 fixed]$ cat pileup.21.0.txt | grep 283422
# 1       283422*a,tq0Q14|:q0Q11
# [wilsonte@gl-login2 fixed]$ cat pileup.21.1.txt | grep 283422
# 1       283422*a,tq0Q11|:q0Q5

# zcat genome_pileup.21.txt.gz
# 21      283422  1       283422*a,t25|:16
# 21      283423  32      :41
# 21      283455  1       :38
# 21      283456  80      :40
# 21      283536  2       283536*tg,aa29|:11
# 21      283538  16      :40
# 21      283554  1       283554*c,g29|:11

# [wilsonte@gl-login2 fixed]$ zcat snv_table.21.txt.gz | grep 283422
# 21      283422  a       t       3       25      25      41

# 281673  281673:39;281712*g,c;281713:28;281741*c,t;281742:32;281774*c,a;281775:78;281853*g,a;281854:110; 111111111       0       1       0       1
# 281673  281673:101;281774*c,a;281775:7;281782*a,g;281783:70;281853*g,a;281854:113;      1110111 0       1       0       1
# 281673  281673:300;     1       0       1       0       13
# 281673  281673:101;281774*c,a;281775:78;281853*g,a;281854:119;  11111   0       1       0       7
# 281673  281673:101;281774*c,a;281775:78;281853*g,a;281854:119;  11101   0       1       0       3
# 281673  281673:76;281749*c,a;281750:223;        101     0       1       0       1
# 281673  281673:101;281774*c,a;281775:7;281782*a,c;281783:70;281853*g,a;281854:113;281967*a,c;281968:5;  111011101       0       1       0       1
# 281673  281673:101;281774*c,a;281775:7;281782*a,c;281783:70;281853*g,a;281854:119;      1110111 0       1       0       1
# 281673  281673:101;281774*c,a;281775:78;281853*g,a;281854:59;281913*c,a;281914:26;281940*c,t;281941:32; 111110101       0       1       0       1
# 281673  281673:39;281712*g,c;281713:28;281741*c,t;281742:32;281774*c,a;281775:2;281777*g,a;281778:75;281853*g,a;281854:119;     11111110111     0       1       0       1
# 281673  281673:39;281712*g,c;281713:28;281741*c,t;281742:32;281774*c,a;281775:78;281853*g,a;281854:119; 111111111       0       1       0       1
# 281673  281673:39;281712*g,c;281713:28;281741*c,t;281742:32;281774*c,a;281775:78;281853*g,a;281854:119; 111111101       0       1       0       1
# 281673  281673:33;281706*t,n;281707:266;        101     0       1       0       1

# 227306  227306:61;227367*,ca;227367:76; 101     0       1       0       1
# 227306  227306:61;227367*,ca;227367:7;227374*t,g;227375:87;     11101   0       1       0       1
# 227306  227306:61;227367*,ca;227367:144;227511*t,c;227512:24;   10111   0       1       0       1
# 227306  227306:61;227367*,ca;227367:144;227511*t,c;227512:92;   10111   0       1       0       1

# 189849  44900817:18;44900815*gg,ttgtt;44900724:91;      111     0       1       0       1
# 189849  44900817:18;44900815*gg,ttgtt;44900710:105;44900709*t,n;44900667:42;    11101   0       1       0       1
# 189849  44900817:18;44900815*gg,ttgtt;44900631:184;     111     0       1       0       1
# 189849  44900817:18;44900815*gg,ttgtt;44900628:187;     111     0       1       0       1

# 837537  44253081:66;44253079*at,gc;44253042:37; 101     0       1       0       1
# 837537  44253081:66;44253079*at,ga;44252865:214;        101     0       1       0       1

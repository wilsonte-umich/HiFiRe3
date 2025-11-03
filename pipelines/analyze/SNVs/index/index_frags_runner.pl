use strict;
use warnings;

# action:
#   fragment indexing for SNV/indel tallying in preparation for pileup

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
require "$ENV{ACTION_DIR}/index/parsing_functions.pl";
require "$ENV{MODULES_DIR}/analyze/SNVs/va_tag_utils.pl";

# environment variables
fillEnvVar(\our $N_CPU,            'N_CPU');
fillEnvVar(\our $SITE_SAM_PREFIX,  'SITE_SAM_PREFIX');
fillEnvVar(\our $GENOME_FASTA,     'GENOME_FASTA');

# shared functions to support SV and SNV/indel indexing
use vars qw($error %chromIndex $snvAlnsH);
our (
    $chrom, $chromIndex1, $paddedChromIndex1, $chromSize,
    @aln, @index,
    $varKey
);

# constants
use constant {
    S_QNAME             => 0, # SITE_SAM fields
    S_FLAG              => 1,
    S_RNAME             => 2,
    S_POS1              => 3,
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    SITE_INDEX1_1       => 6,
    SITE_POS1_1         => 7,
    SITE_DIST_1         => 8,
    SITE_INDEX1_2       => 9,
    SITE_POS1_2         => 10,
    SITE_DIST_2         => 11,
    SEQ_SITE_INDEX1_2   => 12,
    SEQ_SITE_POS1_2     => 13,
    IS_END_TO_END       => 14,
    CH_TAG              => 15,
    TL_TAG              => 16,
    INSERT_SIZE         => 17,
    IS_ALLOWED_SIZE     => 18,
    DE_TAG              => 19,
    HV_TAG              => 20,
    N_REF_BASES         => 21,
    N_READ_BASES        => 22,
    STEM5_LENGTH        => 23,
    STEM3_LENGTH        => 24,
    PASSED_STEM5        => 25,
    PASSED_STEM3        => 26,
    BLOCK_N             => 27,
    ALN_FAILURE_FLAG    => 28,
    JXN_FAILURE_FLAG    => 29,
    TARGET_CLASS        => 30,
    END5_ON_TARGET      => 31,
    READ_HAS_JXN        => 32,
    S_SEQ               => 33,
    S_QUAL              => 34,
    CS_TAG              => 35,
    #-------------
    _IS_PAIRED      => 1,  # SAM FLAG bits
    _PROPER_PAIR    => 2,
    _UNMAPPED       => 4,
    _MATE_UNMAPPED  => 8,
    _REVERSE        => 16,
    _MATE_REVERSE   => 32,
    _FIRST_IN_PAIR  => 64,
    _SECOND_IN_PAIR => 128,
    _SECONDARY      => 256,
    _FAILED_QC      => 512,
    _DUPLICATE      => 1024,
    _SUPPLEMENTAL   => 2048,
    #-------------
    TOP_STRAND    => 0,
    BOTTOM_STRAND => 1,
    #-------------
    INDEX_BIN_SIZE => 1000,
    N_CHAR_ALN_SORT => 18, # leftPos1rightPos1 = 9 + 9
    #-------------
    ALN_HAS_SNV         => 1, # variant flag bits
    READ_HAS_SNV        => 2,
    READ_HAS_SV         => 4,
    READ_HAS_VARIANT    => 8,
};

# initialize the genome
print STDERR "collecting chrom data\n";
setCanonicalChroms();
our @nuclearChroms = getNuclearChroms();
our %chromSizes = getChromSizes("$GENOME_FASTA.fai");

# parallelize by chromosome for speed
# user responsible for setting N_CPU to prevent memory overruns
print STDERR "indexing read alignments by chromosome\n";
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

# parse input SAM
sub processChroms {
    my ($childN) = @_;
    my $readH = $readH[$childN];

    # process data by chrom
    while(my $line = <$readH>){
        chomp $line;

        # open the chromosome-specific site_sam file created by `analyze fragments`
        ($chrom, $chromIndex1, $chromSize) = split("\t", $line);
        my $chromFile = "$SITE_SAM_PREFIX.$chrom.site_sam.gz";

        #############################
        # open my $chromH, "-|", "zcat $chromFile | head -n 100000" or die "$error: could not open: $chromFile: $!\n";
        open my $chromH, "-|", "zcat $chromFile" or die "$error: could not open: $chromFile: $!\n";

        # open caller-specific, chromosome+strand-level output file handles
        # use zero-padded chrom indices to allow simple lexical sort of file names
        $paddedChromIndex1 = sprintf("%02d", $chromIndex1);
        openFileHandles();
        
        # index reads one at a time, counting unique encountered variant patterns per RE fragment
        # indexing action is defined by the calling script
        print STDERR "  $chrom\tassembling index\n";
        @index = ();
        while(my $aln = <$chromH>){
            chomp $aln;
            @aln = split("\t", $aln);
            # filter to non-SV reads (here, based on original alignments in hv tag)
            # thus, there will only ever be one alignment per processed read
            ($aln[IS_END_TO_END] and !($aln[HV_TAG] & READ_HAS_SV)) or next;

            # TODO: do any other work to reverse complement reverse strand alignments
            # everything should leave here oriented as top strand alignments
            # this is fine, since PacBio reads are implicitly duplex and thus not strand-specific

            # parse incoming alignments to stranded VA tag format
            cs_to_va_tag_aln(\@aln);
            # ($aln[S_FLAG] & _REVERSE) and invertReverseStrandAlignment(\@aln, $chromSize); # this is not what we need for HiFiRe3

            # add read to the index
            $index[ int($aln[S_POS1] / INDEX_BIN_SIZE) ]{ indexRead() }++;
        }
        close $chromH;

        # unpack unique encountered variants
        # print to coorindate-sorted output as performed by the calling script
        print STDERR "  $chrom\tsorting and printing bins\n";
        foreach my $binI(0..$#index){
            my $vars = $index[$binI] or next;
            foreach $varKey(sort { 
                substr($a, 0, N_CHAR_ALN_SORT) cmp substr($b, 0, N_CHAR_ALN_SORT) or  # zero-padded leftPos1.rightPos1 sort
                $$vars{$b} <=> $$vars{$a} # most frequent variant pattern sorts first for a given leftPos1.rightPos1
            } keys %$vars){
                printIndexedVariant($$vars{$varKey});
            }
        }
        closeFileHandles();
    }
    reportFinalMetadata_all_thread_chroms($childN);
}

1;

use strict;
use warnings;

# action:
#   fragment indexing for SNV/indel tallying in preparation for pileup
#   indexing of the variant types is handled by the calling scripts

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
require "$ENV{ACTION_DIR}/index/parsing_functions.pl";
require "$ENV{MODULES_DIR}/analyze/SNVs/va_tag_utils.pl";

# environment variables
fillEnvVar(\our $N_CPU,           'N_CPU');
fillEnvVar(\our $SITE_SAM_PREFIX, 'SITE_SAM_PREFIX');
fillEnvVar(\our $GENOME_FASTA,    'GENOME_FASTA');

# shared functions to support SV and SNV/indel indexing
use vars qw($error %chromIndex $snvAlnsH0 $snvAlnsH1);
our (
    $chrom, $chromIndex1, $paddedChromIndex1, $chromSize,
    $prevRead, @alns, @index,
    $strandIndex0, $varKey, $snvAlnsH
);

# constants
use constant {
    S_QNAME             => 0, # SITE_SAM fields
    S_FLAG              => 1,
    S_RNAME             => 2,
    S_POS1              => 3, # 1-based
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    CH_TAG              => 6,
    TL_TAG              => 7,
    DE_TAG              => 8,
    HV_TAG              => 9,
    N_REF_BASES         => 10,
    N_READ_BASES        => 11,
    BLOCK_N             => 12,
    SITE_INDEX1_1       => 13,
    SITE_POS1_1         => 14,
    SITE_DIST_1         => 15,
    SITE_INDEX1_2       => 16,
    SITE_POS1_2         => 17,
    SITE_DIST_2         => 18,
    SEQ_SITE_INDEX1_2   => 19,
    SEQ_SITE_POS1_2     => 20,
    IS_END_TO_END       => 21,
    READ_HAS_JXN        => 22,
    TARGET_CLASS        => 23,
    S_SEQ               => 24,
    S_QUAL              => 25,
    CS_TAG              => 26,
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
        ($prevRead, @alns, @index) = ();
        while(my $aln = <$chromH>){
            chomp $aln;
            my @aln = split("\t", $aln);

            # only process the 5'-most alignment of a read; drop distal, 3' SV alignments
            # S_QNAME has appended readN at this point, so paired reads will always both be processed
            if($prevRead and $prevRead ne $aln[S_QNAME]){
                $strandIndex0 = ($alns[0][S_FLAG] & _REVERSE) ? BOTTOM_STRAND : TOP_STRAND;
                $index[ $strandIndex0 ][ int($alns[0][S_POS1] / INDEX_BIN_SIZE) ]{ indexRead() }++;
                @alns = (); # don't reset @index of course
            }

            # parse incoming alignments to stranded VA tag format and save
            cs_to_va_tag_aln(\@aln);
            ($aln[S_FLAG] & _REVERSE) and invertReverseStrandAlignment(\@aln, $chromSize);
            push @alns, \@aln;
            $prevRead = $aln[S_QNAME];
        }

        # process the last read
        if(@alns){ # `if` here is a catch for siteless chromosomes
            $strandIndex0 = ($alns[0][S_FLAG] & _REVERSE) ? BOTTOM_STRAND : TOP_STRAND;
            $index[ $strandIndex0 ][ int($alns[0][S_POS1] / INDEX_BIN_SIZE) ]{ indexRead() }++;
        }
        close $chromH;

        # unpack unique encountered variants
        # print to coorindate-sorted output as performed by the calling script
        print STDERR "  $chrom\tsorting and printing bins\n";
        foreach $strandIndex0(TOP_STRAND, BOTTOM_STRAND){
            $index[$strandIndex0] or next;
            $snvAlnsH = $strandIndex0 == TOP_STRAND ? $snvAlnsH0 : $snvAlnsH1;
            foreach my $binI(0..$#{$index[$strandIndex0]}){
                my $vars = $index[$strandIndex0][$binI] or next;
                foreach $varKey(sort { 
                    substr($a, 0, N_CHAR_ALN_SORT) cmp substr($b, 0, N_CHAR_ALN_SORT) or  # zero-padded leftPos1.rightPos1 sort
                    $$vars{$b} <=> $$vars{$a} # most frequent variant pattern sorts first for a given leftPos1.rightPos1
                } keys %$vars){
                    printIndexedVariant($$vars{$varKey});
                }
            }
        }
        closeFileHandles();
    }
    reportFinalMetadata_all_thread_chroms($childN);
}

1;

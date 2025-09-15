use strict;
use warnings;

# scope:
#   this script applies to RE-cleaved fragment libraries
# action:
#   shared actions for fragment indexing for SV and SNV/indel tallying
#   indexing of the variant types is handling differently by functions in the calling scripts

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general faidx);
require "$ENV{ACTION_DIR}/index/parsing_functions.pl";

# environment variables
fillEnvVar(\our $N_CPU,                'N_CPU');
fillEnvVar(\our $SITE_SAM_PREFIX,      'SITE_SAM_PREFIX');
fillEnvVar(\our $FILTERING_SITES_FILE, 'FILTERING_SITES_FILE');
fillEnvVar(\our $GENOME_FASTA_SHM,     'GENOME_FASTA_SHM');

# shared functions to support SV and SNV/indel indexing
use vars qw($error %chromIndex @canonicalChroms);
our (
    $chrom, $chromIndex1, $paddedChromIndex1, @sites,
    $prevRead, @alns, @index,
    $siteIndex1, $strandIndex0, $varKey, $faH
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
    SPLIT_TO_S_QUAL     => 27,
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
};

# initialize the genome
print STDERR "collecting chrom index\n";
setCanonicalChroms();
our @nuclearChroms = getNuclearChroms();
loadFaidx($GENOME_FASTA_SHM);

# initialize the conversion from site indices to positions
# once the index is used, it is more informative to track site reference coordinates
print STDERR "loading site positions\n";
open my $sitesH, "-|", "zcat $FILTERING_SITES_FILE" or die "could not open sites file: $FILTERING_SITES_FILE: $!\n";
my $discardHeader = <$sitesH>;
my ($siteIndex1_, $prevChrom) = (0);
while (my $site = <$sitesH>){ # chrom   sitePos1   inSilico   nObserved
    my ($chrom, $sitePos1) = split("\t", $site);
    $prevChrom and $prevChrom ne $chrom and $siteIndex1_ = 0;
    $siteIndex1_++;
    $sites[$chromIndex{$chrom}][$siteIndex1_] = $sitePos1;
    $prevChrom = $chrom;
}
close $sitesH;

# parallelize by chromosome for speed
# user responsible for setting N_CPU to prevent memory overruns
print STDERR "indexing alignments by chromosome\n";
$| = 1;
launchChildThreads(\&processChroms);
use vars qw(@readH @writeH);
my $writeH = $writeH[1];

##############
# foreach my $chrom("chr21", "chr22"){
# foreach my $chrom("chr1"){
foreach my $chrom(@nuclearChroms){
    my $chromIndex1 = $chromIndex{$chrom} or next;
    $writeH = $writeH[$chromIndex1 % $N_CPU + 1];
    print $writeH join("\t", $chrom, $chromIndex1), "\n";
}
finishChildThreads();

# parse input SAM
sub processChroms {
    my ($childN) = @_;
    my $readH = $readH[$childN];

    # open this thread's handle to genome FASTA for pseudo reference contig assembly
    open $faH, "<", $GENOME_FASTA_SHM or die "could not open $GENOME_FASTA_SHM: $!\n";

    # process data by chrom
    while(my $line = <$readH>){
        chomp $line;

        # open the chromosome-specific site_sam file created by `analyze fragments`
        ($chrom, $chromIndex1) = split("\t", $line);
        my $chromFile = "$SITE_SAM_PREFIX.$chrom.site_sam.gz";

        #############################
        # open my $chromH, "-|", "zcat $chromFile | head -n 1000000" or die "$error: could not open: $chromFile: $!\n";
        open my $chromH, "-|", "zcat $chromFile" or die "$error: could not open: $chromFile: $!\n";

        # open caller-specific, chromosome-level output file handles
        # do now in case caller wants to print during indexing
        # use zero-padded chrom indices to allow simple lexical sort of file names
        $paddedChromIndex1 = sprintf("%02d", $chromIndex1);
        openFileHandles();
        
        # index reads one at a time, counting unique encountered variant patterns per RE fragment
        # indexing action is defined by the calling script
        print STDERR "  $chrom\tassembling index\n";
        ($prevRead, @alns, @index) = ();
        while(my $aln = <$chromH>){
            chomp $aln;
            my @aln = split("\t", $aln, SPLIT_TO_S_QUAL);
            if($prevRead and $prevRead ne $aln[S_QNAME]){
                my $key = indexRead();
                $key and $index[ abs($alns[0][SITE_INDEX1_1]) ][ ($alns[0][S_FLAG] & _REVERSE) ? BOTTOM_STRAND : TOP_STRAND ]{ $key }++;
                @alns = (); # don't reset @index of course
            }
            push @alns, \@aln;
            $prevRead = $aln[S_QNAME];
        }
        if(@alns){ # catch for siteless chromosomes
            my $key = indexRead();
            $key and $index[ abs($alns[0][SITE_INDEX1_1]) ][ ($alns[0][S_FLAG] & _REVERSE) ? BOTTOM_STRAND : TOP_STRAND ]{ $key }++;
        }
        close $chromH;

        # unpack unique encountered variants
        # print to output formats as performed by the calling script
        print STDERR "  $chrom\tsorting and printing\n";
        my $i = 1;
        foreach $siteIndex1(1..$#index){
            $index[$siteIndex1] or next;
            foreach $strandIndex0(BOTTOM_STRAND, TOP_STRAND){
                my $vars = $index[$siteIndex1][$strandIndex0] or next;
                foreach $varKey(sort { $$vars{$b} <=> $$vars{$a} } keys %$vars){ # most frequent variant pattern sorts first
                    printIndexedVariant($i, $$vars{$varKey});
                    $i++;
                }
            }
        }
        closeFileHandles();
    }
    reportFinalMetadata_all_thread_chroms($childN);
    close $faH;
}

1;

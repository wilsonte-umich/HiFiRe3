use strict;
use warnings;

# scope:
#   this script is applied to non-RE-cleaved fragments, i.e., a non-ligFree library,
#   to skip non-applicable actions in project_incomplete_reads.pl
# action:
#   record whether:
#       reads have a junction
#       each alignment and read overlapped a target region
# input:
#   name-sorted SITE_SAM on STDIN, with read-level QNAME
# output: 
#   name-sorted SITE_SAM, with read-level QNAME, dispatched to chromosome level files defined by 5'-most alignment
#   where this script:
#       sets:
#           READ_HAS_JXN, which, unlike xf:i: tag, reflects the final read alignment count WITHOUT haplotype override
#           TARGET_CLASS
#   chromosome level output files begin the sorting and indexing process for downstream variant analysis actions

# initialize reporting
our $script = "finish_non_RE_reads";
our $error  = "$script error";
my (@sizeCounts);

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms targets);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\our $SITE_SAM_PREFIX,        'SITE_SAM_PREFIX');
fillEnvVar(\our $INSERT_SIZES_FILE,      'INSERT_SIZES_FILE');
fillEnvVar(\our $TARGETS_BED,            'TARGETS_BED');
fillEnvVar(\our $REGION_PADDING,         'REGION_PADDING');
fillEnvVar(\our $GENOME_SIZE,            'GENOME_SIZE');
my $MAX_INSERT_SIZE = getMaxInsertSize();
our $TARGET_SCALAR = 10;

# initialize the genome and output files
use vars qw(%chromIndex @canonicalChroms);
setCanonicalChroms();

# initialize the target regions, if any
use vars qw($nRegions);
loadTargetRegions();

# constants
use constant {
    S_QNAME             => 0, # SITE_SAM fields
    S_FLAG              => 1,
    S_RNAME             => 2,
    S_POS1              => 3, # 1-based
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    DE_TAG              => 6,
    CS_TAG              => 7,
    XF_TAG              => 8,
    XH_TAG              => 9,
    N_REF_BASES         => 10,
    N_READ_BASES        => 11,
    BLOCK_N             => 12,
    SITE_INDEX1_1       => 13,
    SITE_POS1_1         => 14,
    SITE_HAPS_1         => 15,
    SITE_DIST_1         => 16,
    SITE_INDEX1_2       => 17,
    SITE_POS1_2         => 18,
    SITE_HAPS_2         => 19,
    SITE_DIST_2         => 20,
    SEQ_SITE_INDEX1_2   => 21,
    SEQ_SITE_POS1_2     => 22,
    SEQ_SITE_HAPS_2     => 23,
    IS_END_TO_END       => 24,
    READ_HAS_JXN        => 25,
    TARGET_CLASS        => 26,
    S_SEQ               => 27,
    S_QUAL              => 28,
    #-------------
    SPLIT_TO_TARGET_CLASS => 28,
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
    FALSE  => 0,
    TRUE   => 1,
    # -------------
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
};

# open chrom-level output file handles
my %chromHs;
foreach my $chrom(getNuclearChroms()){
    my $chromFile = "$SITE_SAM_PREFIX.$chrom.site_sam.gz";
    open my $chromH, "|-", "gzip -c > $chromFile" or die "could not open file: $chromFile: $!\n";
    $chromHs{$chrom} = $chromH;
}

# process alignments one read at a time
my ($prevRead, @alns);
while(my $aln = <STDIN>){
    my @aln = split("\t", $aln, SPLIT_TO_TARGET_CLASS);
    if($prevRead and $prevRead ne $aln[S_QNAME]){
        processRead();
        @alns = ();
    }
    push @alns, \@aln;
    $prevRead = $aln[S_QNAME];
}
processRead();

# close chrom output handles
foreach my $chrom(keys %chromHs){
    my $chromH = $chromHs{$chrom};
    close $chromH;
}

# print and save insert size distribution
open my $sizesH, "|-", "gzip -c > $INSERT_SIZES_FILE" or die "could not open pipe: $INSERT_SIZES_FILE: $!\n";
foreach my $size(1..$MAX_INSERT_SIZE){
    print $sizesH join("\t", $size, $sizeCounts[$size] || 0), "\n";
}
close $sizesH;

# report target counts
$nRegions and printTargetCounts($GENOME_SIZE);

# process sequences by QNAME
sub processRead {

    # determine if the read starts on/indexes to a chromosome we intend to analyze downstream
    my $aln5 = $alns[0];
    my $chromH = $chromHs{$$aln5[S_RNAME]};
    $chromH or return;

    # set the alignment target region metadata, i.e., TARGET_CLASS
    $nRegions and setAlnTargetClasses(@alns);

    # update and commit all alignments
    foreach my $i(0..$#alns){
        $alns[$i][READ_HAS_JXN] = @alns > 1 ? TRUE : FALSE;
        print $chromH join("\t", @{$alns[$i]});
    }

    # record the insert size of 5' index alignments for assessing insert size distributions
    # only record read1 since paired read2 carries redundant information
    $$aln5[READ_HAS_JXN] and return;
    my $readN = ($$aln5[S_FLAG] & _IS_PAIRED and $$aln5[S_FLAG] & _SECOND_IN_PAIR) ? READ2 : READ1;
    $readN == READ1 or return;
    # my $size = $$aln5[IS_END_TO_END] ? 
    #     abs($$aln5[SEQ_SITE_POS1_2] - $$aln5[SITE_POS1_1]) + 1 : # for tagFreeLR, these are site pos and not useful for sizing
    #     $$aln5[N_READ_BASES];
    my $size = $$aln5[N_READ_BASES];
    $size > $MAX_INSERT_SIZE and return;
    $sizeCounts[$size]++;
}

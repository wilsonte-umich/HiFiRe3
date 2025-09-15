use strict;
use warnings;

# action:
#   examine 5' most and 3'-most read alignments for a prior site projection of their 3' ends
#   if missing, project them now to the next RE sites
#   record whether:
#       these projections are in a read with a junction, i.e., the nature of the 3' nodes being projected
#       each alignment and read overlapped a target region
#   alignment projection is indicated to fill out RE fragments for plotting and insert size tallies when:
#       a sequence was a truncated, i.e., was an incomplete read on a platform with variable-length single reads
#       a sequence had an SV that disconnected the 5'-most and 3'-most alignments of the parent sequence
#   end projections are used downstream during fragment indexing and path contruction
#   where all outermost endpoints of all sequences are described as the best matching RE site and orientation
# input:
#   name-sorted SITE_SAM on STDIN
# output: 
#   name-sorted SITE_SAM, dispatched to chromosome level files defined by 5'-most alignment
#   where this script:
#       sets:
#           READ_HAS_JXN, which, unlike hv:i: tag, reflects the final read alignment count
#           TARGET_CLASS
#           unmodified alignments either:
#               already had values for SEQ_SITE_INDEX1_2 to SEQ_SITE_POS1_2
#               were internal alignments where projection values stay as default 0
#   chromosome level output files begin the sorting and indexing process for downstream variant analysis actions

# initialize reporting
our $script = "project_incomplete_reads";
our $error  = "$script error";
my @insertSizeCounts;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms targets);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\our $SITE_CHROM_DATA,           'SITE_CHROM_DATA');      # in order of usage: first access the chrom's data
fillEnvVar(\our $SITE_DATA_LOOKUP_WRK,      'SITE_DATA_LOOKUP_WRK'); # then acquire the position and matching haplotypes
fillEnvVar(\our $SITE_SAM_PREFIX,           'SITE_SAM_PREFIX');
fillEnvVar(\our $FILTERED_INSERT_SIZES_FILE,'FILTERED_INSERT_SIZES_FILE');
fillEnvVar(\our $TARGETS_BED,               'TARGETS_BED');
fillEnvVar(\our $REGION_PADDING,            'REGION_PADDING');
fillEnvVar(\our $GENOME_SIZE,               'GENOME_SIZE');
our $TARGET_SCALAR = 10;

# initialize the genome and output files
use vars qw(%chromIndex @canonicalChroms);
setCanonicalChroms();

# initialize the target regions, if any
use vars qw($nRegions);
loadTargetRegions();

# constants
use constant {
    chrom_              => 0, # site index table columns
    chromIndex_         => 1,
    nSites_             => 2,
    chromSize_          => 3,
    closestSiteOffset_  => 4, # one lookup offset per chrom
    siteDataOffset_     => 5,
    #-------------
    SEEK_SET               => 0,   # see perl seek documentation
    CLOSEST_SITE_PACKING   => "l", # signed integer (siteIndex1)
    BYTES_PER_CLOSEST_SITE => 4,
    SITE_DATA_PACKING      => "LSS", # unsigned integer (sitePos1), unsigned shorts (two flanking frag sizes)
    BYTES_PER_SITE_DATA    => 8,
    #-------------
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
    SPLIT_TO_S_SEQ => 26,
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
    P_SITE_INDEX1  => 0,
    P_SITE_POS1    => 1,
    #-------------
    DIR_FORWARD  =>  1,
    DIR_REVERSE  => -1,
    #-------------
    FALSE  => 0,
    TRUE   => 1,
    #-------------
    SIZE_PLOT_BIN_SIZE => 250,
    READ_LEN => 0,
    PROJ_LEN => 1
};

# do not need to load RE site metadata
# we only project sites here by incrementing or decrementing from closest site indices determined upstream

# load the RE site lookup index
my (%chromData);
open my $inH, "<", $SITE_CHROM_DATA  or die "could not open file $SITE_CHROM_DATA : $!";
my $header = <$inH>; # chrom,chromIndex,nSites,chromSize,closestSiteOffset,siteDataOffset
while (my $line = <$inH>){
    chomp $line;
    my @chrom = split("\t", $line);
    $chromData{$chrom[chrom_]} = \@chrom;
}
close $inH;

# open handles to the RE filtering binary site lookups
open my $siteDataH, '<:', $SITE_DATA_LOOKUP_WRK    or die "could not open file: $!";
my $siteDataRaw = "";

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
    my @aln = split("\t", $aln, SPLIT_TO_S_SEQ);
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
open my $sizesH, '>', "$FILTERED_INSERT_SIZES_FILE" or throwError("could not open insert sizes file: $!");
foreach my $sizeType(READ_LEN, PROJ_LEN){
    my $sizeTypeName = $sizeType == READ_LEN ? "filtered" : "projected";
    my @insertSizes = sort {$a <=> $b} keys %{$insertSizeCounts[$sizeType]};
    print $sizesH join("\n", map { 
        join("\t", 
            $sizeTypeName,
            $_, 
            $insertSizeCounts[$sizeType]{$_}
        )
    } @insertSizes), "\n";
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

    # determine if the 5'-most alignment still needs to be projected
    $$aln5[SEQ_SITE_INDEX1_2] or @$aln5[SEQ_SITE_INDEX1_2..SEQ_SITE_POS1_2] = getProjection($aln5);

    # determine if the 3'-most alignment still needs to be projected
    # here projection is more limited and simply bumps to the next RE site on any haplotype
    my $aln3 = $alns[-1];
    $$aln3[SEQ_SITE_INDEX1_2] or @$aln3[SEQ_SITE_INDEX1_2..SEQ_SITE_POS1_2] = getProjection($aln3);

    # set the alignment target region metadata, i.e., TARGET_CLASS
    $nRegions and setAlnTargetClasses(@alns);

    # update and commit all alignments with new or old projection 
    foreach my $i(0..$#alns){
        $alns[$i][READ_HAS_JXN] = @alns > 1 ? TRUE : FALSE;
        print $chromH join("\t", @{$alns[$i]});
    }

    # record the insert size of proper 5' index fragments for assessing insert size distributions
    # goal is to determine the impact of fragment filtering and site projection on insert size distributions
    # as compared to the original unfiltered read lengths obtained prior to alignment
    $$aln5[READ_HAS_JXN] and return;
    my $readLen = length($$aln5[S_SEQ]);
    my $projLen = abs($$aln5[SEQ_SITE_POS1_2] - $$aln5[SITE_POS1_1]) + 1;
    $insertSizeCounts[READ_LEN]{int($readLen / SIZE_PLOT_BIN_SIZE + 0.5) * SIZE_PLOT_BIN_SIZE}++;
    $insertSizeCounts[PROJ_LEN]{int($projLen / SIZE_PLOT_BIN_SIZE + 0.5) * SIZE_PLOT_BIN_SIZE}++;
}
sub getProjection {
    my ($aln) = @_;
    my @proj = @{$aln}[SITE_INDEX1_2..SITE_POS1_2];
    if($$aln[S_FLAG] & _REVERSE){     # no sign on either, pretend an exact match
        $proj[P_SITE_INDEX1] < 0 and  # actual and projected site indices could be the same
        abs($proj[P_SITE_INDEX1]) > 1 and
        $proj[P_SITE_INDEX1] = abs($proj[P_SITE_INDEX1]) - 1;
    } else {
        $proj[P_SITE_INDEX1] > 0 and 
        $proj[P_SITE_INDEX1] < ${$chromData{$$aln[S_RNAME]}}[nSites_] and 
        $proj[P_SITE_INDEX1]++; 
    }
    $proj[P_SITE_INDEX1] < 0 and $proj[P_SITE_INDEX1] = abs($proj[P_SITE_INDEX1]);
    $proj[P_SITE_INDEX1] == abs($$aln[SITE_INDEX1_2]) or
        $proj[P_SITE_POS1] = getSiteData($$aln[S_RNAME], $proj[P_SITE_INDEX1]);
    @proj;
}
sub getSiteData {
    my ($chrom, $siteIndex1) = @_;
    my $siteDataOffset = ${$chromData{$chrom}}[siteDataOffset_] + ($siteIndex1 - 1) * BYTES_PER_SITE_DATA;
    seek($siteDataH, $siteDataOffset, SEEK_SET) or die "could not seek in siteData lookup: $!\n";
    if(!read($siteDataH, $siteDataRaw, BYTES_PER_SITE_DATA)){
        # print STDERR "$error: could not read from siteIndex lookup: $chrom, $adjPos1: $!\n";
        return 0;
    }
    my ($sitePos1) = unpack(SITE_DATA_PACKING, $siteDataRaw);
    $sitePos1; # discard fragment sizes here
}

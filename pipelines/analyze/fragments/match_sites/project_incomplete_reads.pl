use strict;
use warnings;

# action:
#   examine 5' most and 3'-most read alignments for a prior site projection of their 3' ends
#   if missing, project them now to the next (haplotype-consistent) RE sites
#   record whether:
#       these projections are in a read with a junction, i.e., the nature of the 3' nodes being projected
#       each alignment and read overlapped a target region
#   alignment projection is indicated to fill out RE fragments for plotting and insert size tallies when:
#       a sequence was a truncated, i.e., was an incomplete read on a platform with variable-length single reads
#       a paired read was orphaned
#       a sequence had an SV that disconnected the 5'-most and 3'-most alignments of the parent sequence
#   end projections are used downstream during fragment indexing and path contruction
#   where all outermost endpoints of all sequences are described as the best matching RE site and orientation
# input:
#   name-sorted SITE_SAM on STDIN, with read-level QNAME
# output: 
#   name-sorted SITE_SAM, with read-level QNAME, dispatched to chromosome level files defined by 5'-most alignment
#   where this script:
#       sets:
#           READ_HAS_JXN, which, unlike xf:i: tag, reflects the final read alignment count WITHOUT haplotype override
#           TARGET_CLASS
#           unmodified alignments either:
#               already had values for SEQ_SITE_INDEX1_2 to SEQ_SITE_HAPS_2
#               were internal alignments where projection values stay as default 0
#   chromosome level output files begin the sorting and indexing process for downstream variant analysis actions

# initialize reporting
our $script = "project_incomplete_reads";
our $error  = "$script error";
my (@sizeCounts);

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
fillEnvVar(\our $INSERT_SIZES_FILE,         'INSERT_SIZES_FILE');
fillEnvVar(\our $TARGETS_BED,               'TARGETS_BED');
fillEnvVar(\our $REGION_PADDING,            'REGION_PADDING');
fillEnvVar(\our $GENOME_SIZE,               'GENOME_SIZE');
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
    chrom_              => 0, # site index table columns
    chromIndex_         => 1,
    nSites_             => 2,
    chromSize_          => 3,
    closestSiteOffset_  => 4, # one lookup offset per chrom
    siteDataOffset_     => 5,
    #-------------
    SEEK_SET               => 0,   # see perl seek documentation
    CLOSEST_SITE_PACKING   => "l", # signed integer (siteIndex1)
    SITE_DATA_PACKING      => "LCSSSSSS", # unsigned integer (sitePos1), unsigned char (haplotypes), unsigned shorts (haplotype frag sizes)
    BYTES_PER_CLOSEST_SITE => 4,
    BYTES_PER_SITE_DATA    => 17,
    #-------------
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
    REFERENCE_BIT  => 1,
    HAPLOTYPE1_BIT => 2,
    HAPLOTYPE2_BIT => 4,
    #-------------
    P_SITE_INDEX1  => 0,
    P_SITE_POS1    => 1,
    P_SITE_HAPS    => 2,
    #-------------
    DIR_FORWARD  =>  1,
    DIR_REVERSE  => -1,
    #-------------
    FALSE  => 0,
    TRUE   => 1,
    # -------------
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
};
my $allBits = REFERENCE_BIT + HAPLOTYPE1_BIT + HAPLOTYPE2_BIT;
my $genotypeBits = HAPLOTYPE1_BIT + HAPLOTYPE2_BIT;
my @matchingHaps = ($allBits);                                  # 0
$matchingHaps[REFERENCE_BIT]                  = $allBits;       # 1
$matchingHaps[HAPLOTYPE1_BIT]                 = HAPLOTYPE1_BIT; # 2
$matchingHaps[HAPLOTYPE1_BIT + REFERENCE_BIT] = HAPLOTYPE1_BIT; # 3
$matchingHaps[HAPLOTYPE2_BIT]                 = HAPLOTYPE2_BIT; # 4
$matchingHaps[HAPLOTYPE2_BIT + REFERENCE_BIT] = HAPLOTYPE2_BIT; # 5
$matchingHaps[$genotypeBits]                  = $genotypeBits;  # 6
$matchingHaps[$genotypeBits  + REFERENCE_BIT] = $genotypeBits;  # 7

# do not need to load RE site metadata or calculate overhang corrections
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

    # determine if the 5'-most alignment still needs to be projected
    if(!$$aln5[SEQ_SITE_INDEX1_2]){

        # perform an initial projection from the 3' end to next site on any haplotype
        my @proj = getInitialProjection($aln5);

        # ensure that the projected site is on a haplotype consistent with the 5' end
        my $direction = ($$aln5[S_FLAG] & _REVERSE) ? DIR_REVERSE : DIR_FORWARD;
        while(
            $proj[P_SITE_INDEX1] > 1 and 
            $proj[P_SITE_INDEX1] < ${$chromData{$$aln5[S_RNAME]}}[nSites_] and 
            !($proj[P_SITE_HAPS] & $matchingHaps[$$aln5[SITE_HAPS_1]])
        ){
            $proj[P_SITE_INDEX1] += $direction;
            @proj[P_SITE_POS1, P_SITE_HAPS] = getSiteData($$aln5[S_RNAME], $proj[P_SITE_INDEX1]);
        }
        @$aln5[SEQ_SITE_INDEX1_2..SEQ_SITE_HAPS_2] = @proj;
    }

    # determine if the 3'-most alignment still needs to be projected
    # here projection is more limited and simply bumps to the next RE site on any haplotype
    my $aln3 = $alns[$#alns];
    $$aln3[SEQ_SITE_INDEX1_2] or @$aln3[SEQ_SITE_INDEX1_2..SEQ_SITE_HAPS_2] = getInitialProjection($aln3);

    # set the alignment target region metadata, i.e., TARGET_CLASS
    $nRegions and setAlnTargetClasses(@alns);

    # update and commit all alignments with new or old projection 
    foreach my $i(0..$#alns){
        $alns[$i][READ_HAS_JXN] = @alns > 1 ? TRUE : FALSE;
        print $chromH join("\t", @{$alns[$i]});
    }

    # record the insert size of proper 5' index fragments for assessing insert size distributions
    # only record read1 since paired read2 carries redundant information
    $$aln5[READ_HAS_JXN] and return;
    my $readN = ($$aln5[S_FLAG] & _IS_PAIRED and $$aln5[S_FLAG] & _SECOND_IN_PAIR) ? READ2 : READ1;
    my $size = abs($$aln5[SEQ_SITE_POS1_2] - $$aln5[SITE_POS1_1]) + 1;
    $size > $MAX_INSERT_SIZE and return;
    $readN == READ1 and $sizeCounts[$size]++;
}
sub getInitialProjection {
    my ($aln) = @_;
    my @proj = @{$aln}[SITE_INDEX1_2..SITE_HAPS_2];
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
        @proj[P_SITE_POS1, P_SITE_HAPS] = getSiteData($$aln[S_RNAME], $proj[P_SITE_INDEX1]);
    @proj;
}
sub getSiteData {
    my ($chrom, $siteIndex1) = @_;
    my $siteDataOffset = ${$chromData{$chrom}}[siteDataOffset_] + ($siteIndex1 - 1) * BYTES_PER_SITE_DATA;
    seek($siteDataH, $siteDataOffset, SEEK_SET) or die "could not seek in siteData lookup: $!\n";
    if(!read($siteDataH, $siteDataRaw, BYTES_PER_SITE_DATA)){
        # print STDERR "$error: could not read from siteIndex lookup: $chrom, $adjPos1: $!\n";
        return (0, 0);
    }
    my ($sitePos1, $siteHaps) = unpack(SITE_DATA_PACKING, $siteDataRaw);
    ($sitePos1, $siteHaps); # discard fragment sizes here
}

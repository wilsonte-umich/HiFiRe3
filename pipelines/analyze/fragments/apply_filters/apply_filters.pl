use strict;
use warnings;

# action:
#   coordinate actions summarized in apply_filters.sh
# input:
#   name-sorted SAM on STDIN
# output: 
#   name-sorted SITE_SAM
#       ordered from read 5' to 3' end
#       dispatched to chromosome level files defined by 5'-most alignment
#           chromosome level files begin sorting and indexing for downstream variant analysis

# initialize reporting
our $script = "apply_filters";
our $error  = "$script error";
my ($nAlns, $nInputEvents, $nOutputReads, $nSvReads, $nJunctions, $nBlocks, $nFailedTraversal) = (0) x 100;
our (%genomeCounts, %discardCounts, %stemLengthCounts, @nJxnsRejected, %readOutcomes);

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms targets);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general smith_waterman);
resetCountFile();

# environment variables
# genome
fillEnvVar(\our $GENOME_SIZE,               'GENOME_SIZE');
fillEnvVar(\our $IS_COMPOSITE_GENOME,       'IS_COMPOSITE_GENOME');
# fillEnvVar(\our $ZOO_FILTER_LIST,   'ZOO_FILTER_LIST');
# platform and library properties
fillEnvVar(\our $SEQUENCING_PLATFORM,        'SEQUENCING_PLATFORM');
fillEnvVar(\our $PLATFORM_MAX_INSERT_SIZE,   'PLATFORM_MAX_INSERT_SIZE');
fillEnvVar(\our $IS_END_TO_END_READ,         'IS_END_TO_END_READ');
fillEnvVar(\our $READ_PAIR_TYPE,             'READ_PAIR_TYPE');
fillEnvVar(\our $READ_LENGTH_TYPE,           'READ_LENGTH_TYPE');
fillEnvVar(\our $MIN_SELECTED_SIZE,          'MIN_SELECTED_SIZE');
fillEnvVar(\our $SELECTED_SIZE_CV,           'SELECTED_SIZE_CV');
fillEnvVar(\our $MIN_ALLOWED_SIZE,           'MIN_ALLOWED_SIZE');
fillEnvVar(\our $HAS_BASE_ACCURACY,          'HAS_BASE_ACCURACY');
fillEnvVar(\our $EXPECTING_ENDPOINT_RE_SITES,'EXPECTING_ENDPOINT_RE_SITES');
fillEnvVar(\our $REJECTING_JUNCTION_RE_SITES,'REJECTING_JUNCTION_RE_SITES');
# RE site matching
# fillEnvVar(\our $ENZYME_NAME,               'ENZYME_NAME');
# fillEnvVar(\our $BLUNT_RE_TABLE,            'BLUNT_RE_TABLE');
fillEnvVar(\our $SITE_CHROM_DATA_FILE,      'SITE_CHROM_DATA_FILE');     # in order of usage: first access the chrom's data
fillEnvVar(\our $CLOSEST_SITE_LOOKUP_WRK,   'CLOSEST_SITE_LOOKUP_WRK');  # then find the closest site on the chrom to query pos1
fillEnvVar(\our $SITE_DATA_LOOKUP_WRK,      'SITE_DATA_LOOKUP_WRK');     # then acquire the position and matching haplotypes
# tolerances and thresholds
fillEnvVar(\our $CLIP_TOLERANCE,            'CLIP_TOLERANCE');      # RE site matching
fillEnvVar(\our $ACCEPT_ENDPOINT_DISTANCE,  'ACCEPT_ENDPOINT_DISTANCE');
fillEnvVar(\our $REJECT_JUNCTION_DISTANCE,  'REJECT_JUNCTION_DISTANCE');
fillEnvVar(\our $MIN_TRAVERSAL_DELTA,       'MIN_TRAVERSAL_DELTA'); # block numbering
fillEnvVar(\our $MIN_MAPQ,                  'MIN_MAPQ');            # alignment quality
fillEnvVar(\our $MAX_DIVERGENCE,            'MAX_DIVERGENCE');
fillEnvVar(\our $MIN_FLANK_LEN,             'MIN_FLANK_LEN');
fillEnvVar(\our $MIN_AVG_BASE_QUAL,         'MIN_AVG_BASE_QUAL');
# chimeric follow-on and internal adapter parameters
fillEnvVar(\our $INSERTION_WINDOW_SIZE,     'INSERTION_WINDOW_SIZE');
fillEnvVar(\our $MIN_INSERTION_WINDOW_QUAL, 'MIN_INSERTION_WINDOW_QUAL');
fillEnvVar(\our $INSERTION_ADAPTER_SEQUENCE,'INSERTION_ADAPTER_SEQUENCE');
fillEnvVar(\our $MIN_ADAPTER_LENGTH,        'MIN_ADAPTER_LENGTH'); # only look for adapters in junction inserts in this size range
fillEnvVar(\our $MAX_ADAPTER_LENGTH,        'MAX_ADAPTER_LENGTH');
fillEnvVar(\our $MIN_ADAPTER_SCORE,         'MIN_ADAPTER_SCORE');  # empirically determined SW score that signals presence of an adapter in a junction insertion
# target regions
fillEnvVar(\our $TARGETS_BED,               'TARGETS_BED');
fillEnvVar(\our $REGION_PADDING,            'REGION_PADDING');
# output files
fillEnvVar(\our $SITE_SAM_PREFIX,           'SITE_SAM_PREFIX');
fillEnvVar(\our $FILTERED_INSERT_SIZES_FILE,'FILTERED_INSERT_SIZES_FILE');
fillEnvVar(\our $FILTERED_STEM_LENGTHS_FILE,'FILTERED_STEM_LENGTHS_FILE');

# load more dependencies
if($REJECTING_JUNCTION_RE_SITES){
    map { require "$ENV{APPLY_FILTERS_DIR}/$_.pl" } qw(
        match_nodes_to_sites
    );
}
map { require "$ENV{APPLY_FILTERS_DIR}/$_.pl" } qw(
    check_traversal_delta
    enforce_quality_filters
    split_chimeric_reads
);

# set derived variables
use vars qw(%chromData %insertSizeCounts @nAlnsRejected @avgBaseQual);
our $isONT = $SEQUENCING_PLATFORM eq "ONT";
our $isEndToEndPlatform = ($IS_END_TO_END_READ eq "TRUE");
our $isPairedReadPlatform = ($READ_PAIR_TYPE eq "paired");
our $isSizeSelected = ($MIN_SELECTED_SIZE or $MIN_ALLOWED_SIZE);
our $minAllowedSize = $MIN_ALLOWED_SIZE ? $MIN_ALLOWED_SIZE : $MIN_SELECTED_SIZE / (1 + $SELECTED_SIZE_CV); # i.e., 1N, zero if not size-selected
our $maxAllowedSize = $minAllowedSize * 2;
our $TARGET_SCALAR = 10;
our $sizePlotBinSize = $READ_LENGTH_TYPE eq "short" ? 10 : 250; # bp

# initialize the target regions, if any
use vars qw($nRegions);
loadTargetRegions();

# initialize the genome and output files
use vars qw(%chromIndex);
setCanonicalChroms();
my @nuclearChroms = getNuclearChroms();
my %nuclearChroms = map { $_ => 1 } @nuclearChroms;
our @targetChroms = $nRegions ? getTargetChroms() : @nuclearChroms;

# constants
use constant {
    chrom_              => 0, # site index table columns
    chromIndex_         => 1,
    nSites_             => 2,
    chromSize_          => 3,
    closestSiteOffset_  => 4, # one lookup offset per chrom
    siteDataOffset_     => 5,
    #-------------
    QNAME   => 0, # SAM fields
    FLAG    => 1,
    RNAME   => 2,
    POS1    => 3, # 1-based
    MAPQ    => 4,
    CIGAR   => 5,
    RNEXT   => 6,
    PNEXT   => 7,
    TLEN    => 8,
    SEQ     => 9,
    QUAL    => 10,
    TAGS    => 11,
    #-------------
    SPLIT_TO_TAGS => 12, # split leaves TAGS unsplit
    #-------------
    _IS_PAIRED      => 1, # SAM FLAG bits
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
    IS_END_TO_END_READ  => 14,
    IS_END_TO_END_INSERT=> 15,
    NODE_5              => 16,
    NODE_3              => 17,
    CH_TAG              => 18,
    TL_TAG              => 19,
    INSERT_SIZE         => 20,
    IS_ALLOWED_SIZE     => 21,
    FM_TAG              => 22,
    DE_TAG              => 23,
    HV_TAG              => 24,
    N_REF_BASES         => 25,
    N_READ_BASES        => 26,
    STEM5_LENGTH        => 27,
    STEM3_LENGTH        => 28,
    PASSED_STEM5        => 29,
    PASSED_STEM3        => 30,
    BLOCK_N             => 31,
    ALN_FAILURE_FLAG    => 32,
    JXN_FAILURE_FLAG    => 33,
    TARGET_CLASS        => 34,
    IS_ON_TARGET        => 35,
    READ_HAS_JXN        => 36,
    S_SEQ               => 37,
    S_QUAL              => 38,
    CS_TAG              => 39,
    SAM_TAGS            => 40, # used and then discarded
    # -------------
    FALSE   => 0,
    TRUE    => 1,
    #-------------
    KEPT_SEQUENCE           => "valid on-target reads passed to fragment analysis",
    NONCANONICAL_CHROM      => "reads had their 5' alignment on a non-canonical chromosome",
    FAILED_SITE_LOOKUP      => "reads failed site lookup",
    TELOMERIC_FRAGMENT      => "reads were telomeric fragments",
    CLIP_TOO_LARGE          => "reads had a clip that was too large",
    DISTANCE_TOO_LARGE      => "reads had an outer site distance that was too large",
    UNTRIMMED_ONT_ADAPTER   => "ONT reads lacked a 5' adapter trim",
    IS_OFF_TARGET           => "reads had 5' (ONT) or all (other platforms) off-target alignment(s)",
    #-------------
    TRIM5 => 0, 
    TRIM3 => 1,
    # -------------
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
    #-------------
    TOP_STRAND    => "", # signs for integer-encoded nodes
    BOTTOM_STRAND => "-",
    # -------------
    ALN_FAIL_NONE            => 0,
    ALN_FAIL_MAPQ            => 1,
    ALN_FAIL_DIVERGENCE      => 2,
    ALN_FAIL_FLANK_LEN       => 4,
    ALN_FAIL_AVG_BASE_QUAL   => 8,
    # -------------
    JXN_FAIL_NONE            => 0,
    JXN_FAIL_TRAVERSAL_DELTA => 1,
    JXN_FAIL_NONCANONICAL    => 2,
    JXN_FAIL_SITE_MATCH      => 4,
    JXN_FAIL_FOLDBACK_INV    => 8,
    JXN_FAIL_ONT_FOLLOW_ON   => 16, 
    JXN_FAIL_HAS_ADAPTER     => 32,
    #-------------
    AVG_BASE_QUAL_BIN_SIZE => 5,
    #-------------
    READ_LEN_PROPER       => "nonSV.actual",
    PROJ_LEN_PROPER       => "nonSV.projected",
    READ_LEN_CHIMERIC     => "SV.chimeric", # cannot assess projected size of structural variant reads
    READ_LEN_NON_CHIMERIC => "SV.passed",   # these include all SV reads, not just translocations
};
my $READ_LEN_PROPER = READ_LEN_PROPER;
my $PROJ_LEN_PROPER = PROJ_LEN_PROPER;
my $hasPassedJxn  = "at least one passed junction";
my $allFailedJxns = "all failed junctions";
my $hasNoJxns     = "no junctions";
my @nullSites = (0, 0, 0, 0, 0, 0);

# # load any available zoo rejections (sequence that aligned better to another species)
# my %zooRejections;
# if(-f $ZOO_FILTER_LIST){
#     open my $zooH, "-|", "zcat $ZOO_FILTER_LIST" or die "could not open $ZOO_FILTER_LIST: $!\n";
#     while (my $qName = <$zooH>){
#         chomp $qName;
#         $zooRejections{$qName}++;
#     }
#     close $zooH;
# }
# my $nZooRejections = scalar(keys %zooRejections);

# open chrom-level output file handles
my %chromHs;
foreach my $chrom(@targetChroms){
    my $chromFile = "$SITE_SAM_PREFIX.$chrom.site_sam.gz";
    open my $chromH, "|-", "gzip -c > $chromFile" or die "could not open file: $chromFile: $!\n";
    $chromHs{$chrom} = $chromH;
}

# parse input SAM alignments into reads for processing
my ($prevQName, @alnsByReadN, @samAlns, @alns);
while(my $sam = <STDIN>){
    chomp $sam;
    my @sam = split("\t", $sam, SPLIT_TO_TAGS);

    # # execute zoo rejection
    # $zooRejections{$sam[QNAME]} and next;
    $nAlns++;

    # handle sequences as a set of alignments
    if($prevQName and $prevQName ne $sam[QNAME]){
        processEvent();
        @alnsByReadN = ();
    }
    my $readN = ($sam[FLAG] & _IS_PAIRED and $sam[FLAG] & _SECOND_IN_PAIR) ? READ2 : READ1;
    push @{$alnsByReadN[$readN]}, \@sam;
    $prevQName = $sam[QNAME];
}
processEvent();

# close output handles
$REJECTING_JUNCTION_RE_SITES and finishMatchSites();
foreach my $chrom(keys %chromHs){
    my $chromH = $chromHs{$chrom};
    close $chromH;
}

# print summary counts
sub printSummaryCounts {
    my $sep = "-------------\n";
    print STDERR $sep;
    # printCount(commify($nZooRejections),'nZooRejections', 'reads rejected by zoo');

    # counts of all input reads and alignments
    printCount(commify($nAlns),        'nAlns',        'input alignments');
    printCount(commify($nInputEvents), 'nInputEvents', 'input events (single-end reads or read pairs)');
    printCount(commify($nOutputReads), 'nOutputReads', 'output reads (after splitting unmerged read pairs)');
    print STDERR $sep;

    # report counts by genome when applicable
    if($IS_COMPOSITE_GENOME){
        foreach my $genome(sort keys %genomeCounts){
            printCount(commify($genomeCounts{$genome} || 0), "5' read alignments mapped to genome $genome", '');
        }
        print STDERR $sep;
    }

    # counts of discarded reads by reason
    foreach my $reason(
        KEPT_SEQUENCE, 
        NONCANONICAL_CHROM, 
        FAILED_SITE_LOOKUP, TELOMERIC_FRAGMENT, CLIP_TOO_LARGE, DISTANCE_TOO_LARGE, # UNTRIMMED_ONT_ADAPTER
        IS_OFF_TARGET
    ){
        printCount(commify($discardCounts{$reason} || 0), $reason, '');
    }
    print STDERR $sep;

    # counts of alignment-level filter rejections
    printCount(commify($nAlnsRejected[ALN_FAIL_NONE]          || 0), "alns passed all alignment-level filters", '');
    printCount(commify($nAlnsRejected[ALN_FAIL_MAPQ]          || 0), "alns were rejected with MAPQ < $MIN_MAPQ", '');
    printCount(commify($nAlnsRejected[ALN_FAIL_DIVERGENCE]    || 0), "alns were rejected with gap-corrected divergence > $MAX_DIVERGENCE", '');
    printCount(commify($nAlnsRejected[ALN_FAIL_FLANK_LEN]     || 0), "alns were rejected with alignment length < $MIN_FLANK_LEN", '');
    printCount(commify($nAlnsRejected[ALN_FAIL_AVG_BASE_QUAL] || 0), "alns were rejected with alignment avg. base QUAL < $MIN_AVG_BASE_QUAL", '');
    # if(!$HAS_BASE_ACCURACY){
    #     print STDERR "SV alignment avgBaseQuals\n";
    #     for (my $avgBaseQual = 0; $avgBaseQual <= 50; $avgBaseQual += AVG_BASE_QUAL_BIN_SIZE){
    #         my $bin = $avgBaseQual / AVG_BASE_QUAL_BIN_SIZE;
    #         print STDERR join("\t", $avgBaseQual, $avgBaseQual[$bin] || 0), "\n";
    #     }
    # }
    print STDERR $sep;

    # counts of the structures of kept reads
    printCount(commify($nSvReads),  'nSvReads',  'analyzed reads with >1 alignment');
    printCount(commify($nJunctions),'nJunctions','junctions in analyzed SV reads');
    printCount(commify($nFailedTraversal), 'nFailedTraversal', 'junctions failed traversal delta');
    printCount(commify($nBlocks),          'nBlocks',          'output blocks');
    print STDERR $sep;

    # counts of junction-level filter rejections
    printCount(commify($nJxnsRejected[JXN_FAIL_NONE]            || 0), "junctions passed all junction-level filters", '');
    printCount(commify($nJxnsRejected[JXN_FAIL_TRAVERSAL_DELTA] || 0), "junctions were rejected with traversal delta < $MIN_TRAVERSAL_DELTA", '');
    printCount(commify($nJxnsRejected[JXN_FAIL_NONCANONICAL]    || 0), "junctions were rejected with alignment to a non-canonical chromosome", '');
    printCount(commify($nJxnsRejected[JXN_FAIL_SITE_MATCH]      || 0), "junctions were rejected due to junction RE site match within $REJECT_JUNCTION_DISTANCE bp", '');
    printCount(commify($nJxnsRejected[JXN_FAIL_FOLDBACK_INV]    || 0), "junctions were rejected as foldback inversions", '');
    printCount(commify($nJxnsRejected[JXN_FAIL_ONT_FOLLOW_ON]   || 0), "junctions were rejected as ONT follow-on artifacts", '');
    printCount(commify($nJxnsRejected[JXN_FAIL_HAS_ADAPTER]     || 0), "junctions were rejected with detected adapter sequence", '');
    print STDERR $sep;

    # counts of read junction outcomes
    foreach my $outcome($hasNoJxns, $hasPassedJxn, $allFailedJxns){
        printCount(commify($readOutcomes{$outcome} || 0), "reads had $outcome", '');
    }
    print STDERR $sep;

    # target counts
    $nRegions and printTargetCounts($GENOME_SIZE);
    print STDERR $sep;

    # insert size counts and related distributions
    printChimeraSummary();
}
printSummaryCounts();

# determine if read is unmerged pair, process accordingly
our ($isUnmergedPair, $eventHasJxn);
sub processEvent {
    $nInputEvents++;
    $isUnmergedPair = $alnsByReadN[READ2] ? TRUE : FALSE;
    if($isUnmergedPair){
        # collect outer ends of unmerged pairs to define their insert-level properties
        my @pairedAlns = (
            undef,
            # order alignments from read 5' to 3'
            [sort { getQueryStart0(@{$a}[FLAG, CIGAR]) <=> getQueryStart0(@{$b}[FLAG, CIGAR]) } @{$alnsByReadN[READ1]}],
            [sort { getQueryStart0(@{$a}[FLAG, CIGAR]) <=> getQueryStart0(@{$b}[FLAG, CIGAR]) } @{$alnsByReadN[READ2]}]
        );
        $eventHasJxn = (@{$alnsByReadN[READ1]} > 1 or @{$alnsByReadN[READ2]} > 1) ? TRUE : FALSE;
        # then process each read separately
        foreach my $readN(READ1, READ2){
            @samAlns = @{$pairedAlns[$readN]};
            # adjust QNAMEs to indicate read number, i.e., to be unique per paired read
            # from here forward, unmerged paired reads are essentially treated as fixed-length single reads
            # junctions in anomalous gaps are unverifiable and not considered for rare SV analysis
            # pass the paired 5' read end for proper outer node assignment
            foreach my $sam(@samAlns){ $$sam[QNAME] .= "/$readN"; }
            $discardCounts{ processOutputRead($readN == READ1 ? $pairedAlns[READ2][0] : $pairedAlns[READ1][0]) }++;
        }
    } else {
        # order alignments from read 5' to 3'
        @samAlns = sort { getQueryStart0(@{$a}[FLAG, CIGAR]) <=> getQueryStart0(@{$b}[FLAG, CIGAR]) } @{$alnsByReadN[READ1]};
        $eventHasJxn = @samAlns > 1 ? TRUE : FALSE;
        $discardCounts{ processOutputRead() }++;
    }
}

# process QNAME alignment groups, return read rejection reason for counting
sub processOutputRead {
    my ($pairedAln5) = @_;
    $nOutputReads++;

    # only consider reads whose 5' alignment is on a canonical chromosome
    $nuclearChroms{$samAlns[0][RNAME]} or return NONCANONICAL_CHROM;

    # count intergenomic vs intragenomic reads when applicable
    if($IS_COMPOSITE_GENOME){
        my ($chrom5, $genome5) = ($samAlns[0][RNAME] =~ m/(.+)_(.+)/); # e.g., chr1_(hs1)
        $genomeCounts{$genome5}++;
    }

    # initalize site SAM with direct mapping of SAM fields
    # temporarily retain SAM tags for metadata extraction
    @alns = map {
        my @aln = @{$_}[S_QNAME..S_CIGAR];
        @aln[S_SEQ..S_QUAL] = @{$_}[SEQ..QUAL];
        $aln[SAM_TAGS] = "\t".$$_[TAGS]; # facilitate regex matching
        \@aln;
    } @samAlns;

    # check the ends of all alignments for RE site matches
    # discard sequences if outer ends do not match within tolerance
    # adjust outer site positions for clips
    # reject sequences if, for a required outer endpoint:
    #   outer clip is more than CLIP_TOLERANCE bp in length
    #   clip-adjusted endpoint position is more than ACCEPT_ENDPOINT_DISTANCE bp away from closest RE site
    #        *       *       *     proper sitePos1 values
    #   ----|5-----3/5-----3|----  as nodes are numbered for alignments
    #   ----|3-----5/3-----5|----
    my ($isEndToEndRead, @trims);
    my $tl = $alns[0][SAM_TAGS] =~ m/\ttl:Z:(\S+)/ ? $1 : "0,0";
    foreach my $i(0..$#alns){
        if($alns[$i][S_FLAG] & _REVERSE){ # outer endpoints are site matched even when EXPECTING_ENDPOINT_RE_SITES is FALSE if REJECTING_JUNCTION_RE_SITES is TRUE
            @{$alns[$i]}[SITE_INDEX1_1..SITE_DIST_2] = $REJECTING_JUNCTION_RE_SITES ? (
                findClosestSite($alns[$i][S_RNAME], getEnd($alns[$i][S_POS1], $alns[$i][S_CIGAR]) + 1),
                findClosestSite($alns[$i][S_RNAME], $alns[$i][S_POS1])
            ) : @nullSites;
        } else {
            @{$alns[$i]}[SITE_INDEX1_1..SITE_DIST_2] = $REJECTING_JUNCTION_RE_SITES ? ( 
                findClosestSite($alns[$i][S_RNAME], $alns[$i][S_POS1]),
                findClosestSite($alns[$i][S_RNAME], getEnd($alns[$i][S_POS1], $alns[$i][S_CIGAR]) + 1)
            ) : @nullSites;
        }
        if($i == 0){ # read 5' end
            if($REJECTING_JUNCTION_RE_SITES){
                $alns[$i][SITE_INDEX1_1] or return FAILED_SITE_LOOKUP;

                # reject telomeric fragments, i.e., that would cross past the first or last RE site on chrom
                if($alns[$i][S_FLAG] & _REVERSE){
                    abs($alns[$i][SITE_INDEX1_1]) >  1 or return TELOMERIC_FRAGMENT;
                } else {
                    abs($alns[$i][SITE_INDEX1_1]) <  ${$chromData{$alns[$i][S_RNAME]}}[nSites_] or return TELOMERIC_FRAGMENT;
                }
            }

            # get ONT adapter trims once per read
            @trims = split(",", $tl);

            # perform initial end-to-end status determination
            $isEndToEndRead = (
                $isEndToEndPlatform or 
                ($isONT and $trims[TRIM3]) or # TODO: also check from 3' adapter trim on Aviti_1x300 and Ultima, i.e., add tl tag to those platforms
                ($isPairedReadPlatform and !$isUnmergedPair)
            ) ? TRUE : FALSE;
            # check 5' query end, always required to match a RE site
            my $clip5 = 
                ($isONT and !$trims[TRIM5]) ?
                0 : # when ONT untrimmed, clip would include the adapter...
                (($alns[$i][S_FLAG] & _REVERSE) ? getRightClip($alns[$i][S_CIGAR]) : getLeftClip($alns[$i][S_CIGAR]));
            $clip5 > $CLIP_TOLERANCE and return CLIP_TOO_LARGE;
            if($EXPECTING_ENDPOINT_RE_SITES){
                my $dist5 = abs( ($alns[$i][S_FLAG] & _REVERSE) ? $alns[$i][SITE_DIST_1] + $clip5 : $alns[$i][SITE_DIST_1] - $clip5 );
                $dist5 > $ACCEPT_ENDPOINT_DISTANCE and return DISTANCE_TOO_LARGE;
            }
        }
        if($i == $#alns){ # read 3' end
            !$REJECTING_JUNCTION_RE_SITES or $alns[$i][SITE_INDEX1_2] or return FAILED_SITE_LOOKUP;

            # aln3 is guaranteed by 'locate' to be between 1 and nSitesChrom

            # 3' query end if an end-to-end sequence
            # 3' ends on incomplete platforms/sequences are not _required_ to match an RE, but they may be found to do so below
            my $clip3 = 
                ($isONT and !$trims[TRIM3]) ?
                0 :
                (($alns[$i][S_FLAG] & _REVERSE) ? getLeftClip($alns[$i][S_CIGAR]) : getRightClip($alns[$i][S_CIGAR]));
            my $dist3 = abs( ($alns[$i][S_FLAG] & _REVERSE) ? $alns[$i][SITE_DIST_2] - $clip3 : $alns[$i][SITE_DIST_2] + $clip3 );
            if($isEndToEndRead){
                $clip3 > $CLIP_TOLERANCE and return CLIP_TOO_LARGE;
                $EXPECTING_ENDPOINT_RE_SITES and $dist3 > $ACCEPT_ENDPOINT_DISTANCE and return DISTANCE_TOO_LARGE;
            }

            # promote incomplete sequences to end-to-end status if they got to the fragment end
            $EXPECTING_ENDPOINT_RE_SITES and $dist3 <= $ACCEPT_ENDPOINT_DISTANCE and !$isUnmergedPair and $isEndToEndRead = TRUE;
        }
    }
    my $inEndToEndInsert = ($isEndToEndRead or $isPairedReadPlatform) ? TRUE : FALSE;

    # collect initial 3' site projections of each alignment
    my $readHasJxn = @alns > 1 ? TRUE : FALSE;
    foreach my $aln(@alns){
        @$aln[SEQ_SITE_INDEX1_2..SEQ_SITE_POS1_2] = ($isEndToEndRead and !$readHasJxn) ?
            @{$alns[-1]}[SITE_INDEX1_2..SITE_POS1_2] :
            (0,0);
        $$aln[SEQ_SITE_INDEX1_2] = abs($$aln[SEQ_SITE_INDEX1_2]);
    }

    # determine if the 5'-most alignment still needs to be projected
    my ($node5, $node3);
    if($EXPECTING_ENDPOINT_RE_SITES){
        $alns[ 0][SEQ_SITE_INDEX1_2] or @{$alns[ 0]}[SEQ_SITE_INDEX1_2..SEQ_SITE_POS1_2] = getProjection($alns[0]);

        # determine if the 3'-most alignment still needs to be projected
        # here projection is more limited and simply bumps to the next RE site on any haplotype
        $alns[-1][SEQ_SITE_INDEX1_2] or @{$alns[-1]}[SEQ_SITE_INDEX1_2..SEQ_SITE_POS1_2] = getProjection($alns[-1]);

        # get outer nodes based on RE site positions
        my $pairedSitePos1;
        if($pairedAln5){
            $pairedSitePos1 = (findClosestSite( # may not have assigned site position to other read yet for first read of unmerged pair
                $$pairedAln5[S_RNAME], 
                ($$pairedAln5[S_FLAG] & _REVERSE) ? getEnd($$pairedAln5[S_POS1], $$pairedAln5[S_CIGAR]) : $$pairedAln5[S_POS1],
            ))[1] or return FAILED_SITE_LOOKUP;
        }
        $node5 = parseSignedNode( # reads always define their own 5' outer node
            $alns[0], 
            $alns[0][SITE_POS1_1], # using site positions enables fuzzy endpoint matching when comparing reads
            BOTTOM_STRAND, TOP_STRAND
        );
        $node3 = $pairedAln5 ? 
            parseSignedNode(
                $pairedAln5, 
                $pairedSitePos1,
                TOP_STRAND, BOTTOM_STRAND # invert the strand orientation of paired 5' ends to match the behavior of single or merged reads
            ) : 
            parseSignedNode(
                $alns[-1], 
                $alns[-1][SEQ_SITE_POS1_2],  # using projected 3' ends allows best deduplication of partially sequenced RE fragments
                BOTTOM_STRAND, TOP_STRAND
            );

    # get outer nodes based on alignment positions
    } else {
        $node5 =  parseSignedNode(
            $alns[0], 
            ($alns[0][S_FLAG] & _REVERSE) ? getEnd($alns[0][S_POS1], $alns[0][S_CIGAR]) : $alns[0][S_POS1], 
            BOTTOM_STRAND, TOP_STRAND
        );
        $node3 = $pairedAln5 ? 
            parseSignedNode(
                $pairedAln5, 
                ($$pairedAln5[S_FLAG] & _REVERSE) ? getEnd($$pairedAln5[S_POS1], $$pairedAln5[S_CIGAR]) : $$pairedAln5[S_POS1],
                TOP_STRAND, BOTTOM_STRAND
            ) : 
            parseSignedNode( # this unprojectable non-RE node position may not be the end of the insert if a single read is not end-to-end
                $alns[-1], 
                ($alns[-1][S_FLAG] & _REVERSE) ? $alns[-1][S_POS1] : getEnd($alns[-1][S_POS1], $alns[-1][S_CIGAR]), 
                BOTTOM_STRAND, TOP_STRAND
            )
    }

    # collect next chunk of read-level metadata for passing into each alignment
    my $ch = $alns[0][SAM_TAGS] =~ m/\tch:i:(\d+)/ ? $1 : 0;
    my $readLen = length($alns[0][S_SEQ]);
    my $insertSize = $isUnmergedPair ? -$node3 - $node5 : $readLen;
    $isUnmergedPair and (
        $eventHasJxn or 
        $insertSize < 0 or
        $insertSize > $PLATFORM_MAX_INSERT_SIZE or 
        $alns[0][S_RNAME] ne $$pairedAln5[S_RNAME]
    ) and $insertSize = 0; # when read carries any type of SV, whether sequenced or not, we cannot infer the true insert size
    my $readStart0 = getQueryStart0($alns[ 0][S_FLAG], $alns[ 0][S_CIGAR]);
    my $readEnd1   = getQueryEnd1(  $alns[-1][S_FLAG], $alns[-1][S_CIGAR], $readLen); # may not be the end of the insert
    my $isAllowedSize = (!$isSizeSelected or ($insertSize >= $minAllowedSize and $insertSize <= $maxAllowedSize)) ? TRUE : FALSE;
    foreach my $aln(@alns){
        @$aln[IS_END_TO_END_READ..IS_ALLOWED_SIZE] = (
            $isEndToEndRead, 
            $inEndToEndInsert,
            $node5, # like other read and insert-level properties, these repeat the same on all alignments
            $node3, 
            $ch, 
            $tl, 
            $insertSize, 
            $isAllowedSize
        );
    }

    # stop processing off-target reads, we have all required counts and they won't be used for variant calling
    my $isOnTarget = $nRegions ? setAlnTargetClasses(@alns) : TRUE;
    my $chromH = $chromHs{$samAlns[0][RNAME]};
    ($isOnTarget and $chromH) or return IS_OFF_TARGET;

    # add alignment-level metadata
    my $blockN = 1;
    foreach my $aln(@alns){
        my $stem5Len = getQueryEnd1($$aln[S_FLAG], $$aln[S_CIGAR], $readLen) - $readStart0;
        my $stem3Len = $isEndToEndRead ? $readEnd1 - getQueryStart0($$aln[S_FLAG], $$aln[S_CIGAR]) : 0;
        my $fm = $$aln[SAM_TAGS] =~ m/\tfm:Z:(\S+)/ ? $1 : "0:0:0";
        my $de = $$aln[SAM_TAGS] =~ m/\tde:f:(\S+)/ ? $1 : 0; # do not move these below
        my $hv = $$aln[SAM_TAGS] =~ m/\thv:i:(\d+)/ ? $1 : 0;
        my $cs = $$aln[SAM_TAGS] =~ m/\tcs:Z:(\S+)/ ? $1 : "*";
        @$aln[
            FM_TAG..READ_HAS_JXN,
            CS_TAG
        ] = (
            $fm,
            $de,
            $hv,
            getRefSpan($$aln[S_CIGAR]),     # N_REF_BASES
            getAlignedSize($$aln[S_CIGAR]), # N_READ_BASES
            $stem5Len,
            $stem3Len,
            (!$isSizeSelected or  $stem5Len < $minAllowedSize)                    ? TRUE : FALSE,
            (!$isSizeSelected or ($stem3Len < $minAllowedSize and $stem3Len > 0)) ? TRUE : FALSE,
            $blockN,
            getAlnFailureFlag($aln, $de, $readHasJxn),
            0, # JXN_FAILURE_FLAG
            $$aln[TARGET_CLASS] || 0,
            $isOnTarget,
            0, # READ_HAS_JXN, NOT the same as $readHasJxn = @alns > 1, this will require JXN_FAILURE_FLAG == 0
            #--------------------
            $cs
        );
    }

    # check junction read traversal delta for SV reads
    my $anyJxnPassed = FALSE;
    if(@alns > 1){
        $nSvReads++;
        $nJunctions += @alns - 1;

        # assess traversal delta for every pair of alignments flanking every possible read sub-path
        # fail all junctions within a failed path
        # once a junction fails, it stays failed even if it passes a different pair of alignments
        my @failedTraversalDelta = (0) x @alns;
        foreach my $aln2I(1..$#alns){
            foreach my $aln1I(0..($aln2I - 1)){
                failedTraversalDelta($alns[$aln1I], $alns[$aln2I], $readLen) and
                    @failedTraversalDelta[($aln1I + 1)..$aln2I] = (1) x ($aln2I - $aln1I);
            }
        }

        # update BLOCK_N based on failedTraversalDelta
        foreach my $aln2I(1..$#alns){
            $nFailedTraversal += $failedTraversalDelta[$aln2I];
            $failedTraversalDelta[$aln2I] or $blockN++; # thus, increment alignment blockN on every non-failed junction
            $alns[$aln2I][BLOCK_N] = $blockN;

            # apply chimeric junction filters
            $alns[$aln2I - 1][JXN_FAILURE_FLAG] = $failedTraversalDelta[$aln2I] ? 
                JXN_FAIL_TRAVERSAL_DELTA : 
                getJxnFailureFlag($alns[$aln2I - 1], $alns[$aln2I], $readLen);
            $nJxnsRejected[$alns[$aln2I - 1][JXN_FAILURE_FLAG]]++;
            $alns[$aln2I - 1][JXN_FAILURE_FLAG] or $anyJxnPassed = TRUE;
        }
    }
    $nBlocks += $blockN;

    # collect the distributions of insert sizes of on-target reads
    my $sizeBin = int($insertSize / $sizePlotBinSize) * $sizePlotBinSize;
    if($readHasJxn){ # disregard multi-junction reads when assessing insert sizes
        @alns == 2 and $sizeBin and recordVariantSizes($sizeBin, @alns); 
        $readOutcomes{$anyJxnPassed ? $hasPassedJxn : $allFailedJxns}++;
    } else {
        $sizeBin and $insertSizeCounts{$READ_LEN_PROPER}{$sizeBin}++;
        my $projLen = $EXPECTING_ENDPOINT_RE_SITES ? 
            abs($alns[0][SEQ_SITE_POS1_2] - $alns[0][SITE_POS1_1]) + 1 :
            $insertSize;
        $projLen and $insertSizeCounts{$PROJ_LEN_PROPER}{int($projLen / $sizePlotBinSize) * $sizePlotBinSize}++;
        $readOutcomes{$hasNoJxns}++;
    }

    # print all alignment and return success
    foreach my $aln(@alns){
        $$aln[READ_HAS_JXN] = $anyJxnPassed;
        print $chromH join("\t", @$aln[S_QNAME..CS_TAG]), "\n";
    }

    return KEPT_SEQUENCE;
}

# parse signed node representation for a read outer end
sub parseSignedNode {
    my ($aln, $pos1, $bottomStrand, $topStrand) = @_;
    my $strandSign = ($$aln[S_FLAG] & _REVERSE) ? $bottomStrand : $topStrand;
    $strandSign.$chromIndex{$$aln[S_RNAME]}.(sprintf "%09d", $pos1);
}

1;

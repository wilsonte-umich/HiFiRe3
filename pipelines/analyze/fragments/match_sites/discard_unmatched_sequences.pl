use strict;
use warnings;

# action:
#   discard entire sequences when either complete outermost node does not match a RE site within ACCEPT_ENDPOINT_DISTANCE
# input:
#   name-sorted SITE_SAM on STDIN
# output: 
#   name-sorted SITE_SAM on STDOUT
#   where this script:
#       updates SEQ_SITE_INDEX1_2 to IS_END_TO_END when known
#       appends readN to QNAME

# initialize reporting
our $script = "discard_unmatched_sequences";
our $error  = "$script error";
my ($nSequences) = (0) x 10;
my %discardCounts;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\our $SEQUENCING_PLATFORM,     'SEQUENCING_PLATFORM');
fillEnvVar(\our $READ_PAIR_TYPE,          'READ_PAIR_TYPE');
fillEnvVar(\our $IS_END_TO_END_READ,      'IS_END_TO_END_READ');
fillEnvVar(\our $CLIP_TOLERANCE,          'CLIP_TOLERANCE');
fillEnvVar(\our $ACCEPT_ENDPOINT_DISTANCE,'ACCEPT_ENDPOINT_DISTANCE');
fillEnvVar(\our $SITE_CHROM_DATA,         'SITE_CHROM_DATA');

# set platform-specific parameters
my $isONT = $SEQUENCING_PLATFORM eq "ONT";
my $isPairedReadPlatform = ($READ_PAIR_TYPE eq "paired");
my $isEndToEndPlatform = ($IS_END_TO_END_READ eq "TRUE");

# constants
use constant {
    chrom_              => 0, # site index table columns
    chromIndex_         => 1,
    nSites_             => 2,
    chromSize_          => 3,
    closestSiteOffset_  => 4, # one lookup offset per chrom
    siteDataOffset_     => 5,
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
    SPLIT_TO_IS_END_TO_END => 26,
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
    channel             => 0, # incoming qName extensions
    trim5               => 1,
    trim3               => 2,
    isMerged            => 3, # true=2 for legacy reasons
    nRead1              => 4,
    nRead2              => 5,
    splitGapReadN       => 6,
    N_QNAME_EXTENSIONS  => 7,
    #-------------
    KEPT_SEQUENCE           => "valid sequences passed to fragment analysis",    
    FAILED_SITE_LOOKUP      => "failed site lookup",
    TELOMERIC_FRAGMENT      => "telomeric fragment",
    CLIP_TOO_LARGE          => "clip too large",
    DISTANCE_TOO_LARGE      => "distance too large",
    UNTRIMMED_ONT_ADAPTER   => "ONT without 5' adapter trim",
    # -------------
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
    # -------------
    FALSE  => 0,
    TRUE   => 1,
};

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

# process alignments one at a time
my (@alns, $prevQName);
while(my $aln = <STDIN>){
    my @aln = split("\t", $aln, SPLIT_TO_IS_END_TO_END);
    if($prevQName and $prevQName ne $aln[S_QNAME]){
        $discardCounts{ processQName() }++;
        @alns = ();
    }
    my $readN = ($aln[S_FLAG] & _IS_PAIRED and $aln[S_FLAG] & _SECOND_IN_PAIR) ? READ2 : READ1;
    push @{$alns[$readN]}, \@aln;
    $prevQName = $aln[S_QNAME];
}
$discardCounts{ processQName() }++;

# print summary information
printCount(commify($nSequences), 'nSequences', 'input sequences');
foreach my $reason(KEPT_SEQUENCE, FAILED_SITE_LOOKUP, TELOMERIC_FRAGMENT, CLIP_TOO_LARGE, DISTANCE_TOO_LARGE){ # UNTRIMMED_ONT_ADAPTER
    printCount(commify($discardCounts{$reason} || 0), $reason, '');
}

# process sequences by QNAME
sub processQName {
    $nSequences++;

    # parse to two outermost alignments of the query sequence
    # 5 and 3 refer to the 5' and 3' ends of the query molecule as sequenced
    # 3' end is found in different ways depending on pairing status
    my $aln5 = $alns[READ1][0];
    my $aln3 = $alns[READ2] ? $alns[READ2][0] : $alns[READ1][$#{$alns[READ1]}];
    $$aln5[SITE_INDEX1_1] or return FAILED_SITE_LOOKUP;
    if($alns[READ2]){
        $$aln3[SITE_INDEX1_1] or return FAILED_SITE_LOOKUP;
    } else {
        $$aln3[SITE_INDEX1_2] or return FAILED_SITE_LOOKUP;
    }
    
    # require that ONT reads have trimmed 5' adapters
    my @qName = split(":", $$aln5[S_QNAME]);
    my @extensions = splice(@qName, -N_QNAME_EXTENSIONS);
    # $isONT and !$extensions[trim5] and return UNTRIMMED_ONT_ADAPTER; DEPRECATED, see below

    # reject telomeric fragments, i.e., that would cross past the first or last RE site on chrom
    if($$aln5[S_FLAG] & _REVERSE){
        abs($$aln5[SITE_INDEX1_1]) >  1 or return TELOMERIC_FRAGMENT;
    } else {
        abs($$aln5[SITE_INDEX1_1]) <  ${$chromData{$$aln5[S_RNAME]}}[nSites_] or return TELOMERIC_FRAGMENT;
    }
    # aln3 is guaranteed by 'locate' to be between 1 and nSitesChrom

    # adjust outer site positions for clips
    # reject sequences if, for a required outer endpoint:
    #   outer clip is more than CLIP_TOLERANCE bp in length
    #   clip-adjusted endpoint position is more than ACCEPT_ENDPOINT_DISTANCE bp away from closest RE site
    #        *       *       *     proper sitePos1 values
    #   ----|5-----3/5-----3|----  as nodes are numbered for alignments
    #   ----|3-----5/3-----5|----

    # 5' query end, always required to match a RE site
    my $clip5 = 
        ($isONT and !$extensions[trim5]) ?
        0 : # when ONT untrimmed, clip would include the adapter...
        (($$aln5[S_FLAG] & _REVERSE) ? getRightClip($$aln5[S_CIGAR]) : getLeftClip($$aln5[S_CIGAR]));
    $clip5 > $CLIP_TOLERANCE and return CLIP_TOO_LARGE;
    my $dist5 = abs( ($$aln5[S_FLAG] & _REVERSE) ? $$aln5[SITE_DIST_1] + $clip5 : $$aln5[SITE_DIST_1] - $clip5 );
    $dist5 > $ACCEPT_ENDPOINT_DISTANCE and return DISTANCE_TOO_LARGE;

    # 3' query end if an end-to-end sequence
    # 3' ends on incomplete platforms/sequences are not _required_ to match an RE, but they may be found to do so below
    my $isEndToEnd = (
        $isEndToEndPlatform 
        or 
        (
            $isPairedReadPlatform and 
            (
                $extensions[isMerged] or
                $alns[READ2]
                # remaining paired are orphaned reads
            )
        ) 
        or 
        (
            $isONT and 
            $extensions[trim3]
            # incomplete single reads remain
        )
    ) ? TRUE : FALSE;
    my ($clip3, $dist3); # the 3' end of the query sequence, so could be the 5' end of paired read 2
    if($alns[READ2]){
        $clip3 = ($$aln3[S_FLAG] & _REVERSE) ? getRightClip($$aln3[S_CIGAR]) : getLeftClip($$aln3[S_CIGAR]);
        $dist3 = abs( ($$aln3[S_FLAG] & _REVERSE) ? $$aln3[SITE_DIST_1] + $clip3 : $$aln3[SITE_DIST_1] - $clip3 );
    } else {
        $clip3 = 
            ($isONT and !$extensions[trim3]) ?
            0 :
            (($$aln3[S_FLAG] & _REVERSE) ? getLeftClip($$aln3[S_CIGAR]) : getRightClip($$aln3[S_CIGAR]));
        $dist3 = abs( ($$aln3[S_FLAG] & _REVERSE) ? $$aln3[SITE_DIST_2] - $clip3 : $$aln3[SITE_DIST_2] + $clip3 );
    }
    if($isEndToEnd){
        $clip3 > $CLIP_TOLERANCE and return CLIP_TOO_LARGE;
        $dist3 > $ACCEPT_ENDPOINT_DISTANCE and return DISTANCE_TOO_LARGE;
    }

    # promote incomplete sequences to end-to-end status if they got to the fragment end
    $dist3 <= $ACCEPT_ENDPOINT_DISTANCE and $isEndToEnd = TRUE;

    # append readN to all alignments and commit to next actions, which all act at read level
    my $hasSV = ( # do not use EVENT_HAS_SV here, need to establish connection between ends according to reference
        @{$alns[READ1]} > 1 or
        (
            $alns[READ2] and 
            @{$alns[READ2]} > 1  # gaps cannot and do not have junctions, as enforced by order_alignments
        )
    );
    foreach my $readN(READ1, READ2){
        $alns[$readN] or next;
        foreach my $aln(@{$alns[$readN]}){
            $$aln[S_QNAME] .= ":$readN";
            if($isEndToEnd and !$hasSV){ 
                @$aln[SEQ_SITE_INDEX1_2..SEQ_SITE_HAPS_2] = 
                    $readN == READ1 ? 
                    (
                        $alns[READ2] ? 
                        @$aln3[SITE_INDEX1_1..SITE_HAPS_1] :
                        @$aln3[SITE_INDEX1_2..SITE_HAPS_2]
                    ) : 
                    @$aln5[SITE_INDEX1_1..SITE_HAPS_1];
            }
            $$aln[SEQ_SITE_INDEX1_2] = abs($$aln[SEQ_SITE_INDEX1_2]);
            $$aln[IS_END_TO_END] = $isEndToEnd;
            print join("\t", @$aln);
        }
    }
    return KEPT_SEQUENCE;
}

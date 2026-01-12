use strict;
use warnings;

# action:
#   find ONT adapters on previously untrimmed read ends
#   search includes blunt RE half-cut-site to enhance search power vs. adapter alone
#   trim just the adapter portion, leaving the half-cut-site genomic bases intact
#   add tag indicating trim status
# input:
#   an unaligned SAM stream, specifically, a Dorado-compatible output stream, on STDIN
# output:
#   updated SAM stream on STDOUT, with new tags tl:Z: = adapter trim lengths in format `<5' trim>,<3' trim>`
# notes:
#   this script processes unaligned ONT reads, i.e., never reverse complemented, so first bases correspond to the 5' adapter
#   all reads are expected to have 5' adapters, but some may be low quality bases, truncated, or chimeras split by Dorado
#   3' adapters are always shorter, often truncated, and not always present for incompletely sequenced fragments

# load dependencies
our $script = "trim_adapters";
our $error  = "$script error";
our ($matchScore, $mismatchPenalty, $gapOpenPenalty, $gapExtensionPenalty) = 
    (1,           -1,               -1.5,            -1); # penalties lowered from defaults of -1.5, -2.5, -1
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general smith_waterman);

# environment variables
fillEnvVar(\our $BLUNT_RE_TABLE,   'BLUNT_RE_TABLE'); # HiFiRe3 only supports blunt REs
fillEnvVar(\our $ENZYME_NAME,      'ENZYME_NAME');
fillEnvVar(\our $ADAPTER_SEQUENCE, 'ADAPTER_SEQUENCE');
fillEnvVar(\our $PLATFORM_MIN_INSERT_SIZE, 'PLATFORM_MIN_INSERT_SIZE');

# set operating parameters
my $MIN_SCORE = 10; # 3 RE bases + 1 A-tail + 6 adapter bases
my $includeREHalfSite = ($ENZYME_NAME ne "NA");

# constants
use constant {
    QNAME   => 0, # SAM fields
    FLAG    => 1,
    RNAME   => 2,
    POS     => 3, # 1-based
    MAPQ    => 4,
    CIGAR   => 5,
    RNEXT   => 6,
    PNEXT   => 7,
    TLEN    => 8,
    SEQ     => 9,
    QUAL    => 10,
    TAGS    => 11,
    #--------------
    NOT_FAST        => 0,
    FORCE_QRY_START => 1,  # smith_waterman control
    FORCE_QRY_END   => 2,
    SEARCH_SPACE_5  => 60, # search spaces empirically determined from actual reads
    SEARCH_SPACE_3  => 20, # wide enough to be sensitive, narrow enough to be fast and specific
};

# initialize the RE
my ($cutSiteLeftHalf, $cutSiteRightHalf, $cutSiteOffset) = ("", "", 0);
if($includeREHalfSite){
    open my $inH, "<", $BLUNT_RE_TABLE or die "could not open: $BLUNT_RE_TABLE: $!\n";
    my $header = <$inH>; # enzyme,strand,cut_site,regex,offset,CpG_priority,...
    while (my $line = <$inH>){
        my ($enzyme, $strand, $cut_site, $regex, $offset) = split(",", $line);
        $enzyme eq $ENZYME_NAME or next;
        $cut_site = uc($cut_site);
        $cutSiteLeftHalf  = substr($cut_site, 0, $offset);
        $cutSiteRightHalf = substr($cut_site, $offset, $offset);
        $cutSiteOffset    = $offset;
    }
    close $inH;
    $cutSiteLeftHalf or die "unrecognized enzyme; must be a blunt cutter in shared/modules/REs/blunt_enzymes.csv\n";
}

# initalize the adapter search sequences
my $adapterCore = $ADAPTER_SEQUENCE; # duplex portion of the ONT kit adapter; for ligation kit, last T matches the one-base A-tail; fused to 5' genomic ends 
my $adapter5 = {
    adapter   => $adapterCore,
    adapterRE => $adapterCore.$cutSiteRightHalf
};
$$adapter5{exact} = substr($$adapter5{adapterRE}, -$MIN_SCORE);
my $adapterCoreRc = $adapterCore; # fused to 3' genomic ends
rc(\$adapterCoreRc);
my $adapter3 = {
    adapter   => $adapterCoreRc,
    adapterRE => $cutSiteLeftHalf.$adapterCoreRc
};
$$adapter3{exact} = substr($$adapter3{adapterRE}, 0, $MIN_SCORE);

# trim reads one at a time, on both ends
while (my $read = <STDIN>){
    if($read =~ m/^\@/){
        print $read;
    } else {
        chomp $read;
        my @read = split("\t", $read, 12);
        length($read[SEQ]) >= $PLATFORM_MIN_INSERT_SIZE or next; # reject reads too small to even search for adapaters
        my ($found5, $trimLen5) = find_adapter($adapter5, substr($read[SEQ], 0, SEARCH_SPACE_5), FORCE_QRY_END);
        my ($found3, $trimLen3) = find_adapter($adapter3, substr($read[SEQ], -SEARCH_SPACE_3),   FORCE_QRY_START);
        if($found3){ # order is important
            $read[SEQ]  = substr($read[SEQ],  0, -$trimLen3);
            $read[QUAL] = substr($read[QUAL], 0, -$trimLen3);
        }
        if($found5){ # order is important
            $read[SEQ]  = substr($read[SEQ],  $trimLen5);
            $read[QUAL] = substr($read[QUAL], $trimLen5);
        }
        # add trim tag and print
        print join("\t", @read, "tl:Z:".($found5 ? $trimLen5 : 0).",".($found3 ? $trimLen3 : 0)), "\n";
    }
}

# execute an efficient adapter search
# searches force alignment to the end of the adapter closest to the adapter-gDNA transition for specificity (also the highest quality bases)
sub find_adapter {
    my ($adapter, $readSegment, $forcedEnd) = @_;
    my ($qryOnRef, $score, $startQry0, $endQry0, $startRef0, $endRef0, $trimLen);

        # for exploring and debugging smith_watermen
        # ($qryOnRef, $score, $startQry0, $endQry0, $startRef0, $endRef0) = 
        #     smith_waterman($$adapter{adapterRE}, $readSegment, NOT_FAST, $forcedEnd);
        # $score or return;
        # my @ref = split("", $readSegment);
        # @ref = map { $ref[$_ + $startRef0].("-" x (length($$qryOnRef[$_]) - 1)) } 0..$#$qryOnRef;
        # print join("\n",
        #     (" " x $startRef0).join("", @$qryOnRef),
        #     substr($readSegment, 0, $startRef0).join("", @ref).substr($readSegment, $endRef0 + 1),
        #     $score
        # ), "\n\n";
        # return;

    # for speed, search first for the adapter-RE transition by exact matching
    $startRef0 = index($readSegment, $$adapter{exact});
    if($startRef0 != -1){
        $score = $MIN_SCORE;
        $trimLen = $forcedEnd == FORCE_QRY_END ? # as measured from SEQ 5' or 3' end; corresponds to the A-tail base
            $startRef0 + $MIN_SCORE - $cutSiteOffset : 
            SEARCH_SPACE_3 - $startRef0 - $cutSiteOffset;

    # if not found, use smith_waterman to locate the adapter sequence with the RE half-site
    # will have a sequencing error within the 10bp exact search region
    } else {
        ($qryOnRef, $score, $startQry0, $endQry0, $startRef0, $endRef0) = 
            smith_waterman($$adapter{adapterRE}, $readSegment, NOT_FAST, $forcedEnd);
        if($score >= $MIN_SCORE){
            $trimLen = $forcedEnd == FORCE_QRY_END ? 
                $endRef0 + 1 - $cutSiteOffset: 
                SEARCH_SPACE_3 - $startRef0 - $cutSiteOffset
        
        # if not found, use smith_waterman to locate an adapter-gDNA transition without the RE half-site
        # in addition to sequencing errors, this may arise from adapter ligation onto a non-RE gDNA end
        } elsif ($includeREHalfSite) {
            ($qryOnRef, $score, $startQry0, $endQry0, $startRef0, $endRef0) = 
                smith_waterman($$adapter{adapter}, $readSegment, NOT_FAST, $forcedEnd);
            if($score >= $MIN_SCORE){
                $trimLen = $forcedEnd == FORCE_QRY_END ?
                    $endRef0 + 1: 
                    SEARCH_SPACE_3 - $startRef0
            }
        }
    }
    return ($score >= $MIN_SCORE, $trimLen);
}

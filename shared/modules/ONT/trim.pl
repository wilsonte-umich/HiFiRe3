use strict;
use warnings;

# action:
#   find ONT adapters on previously untrimmed read ends
#   for ligFree, search includes RE half-cut-site to enhance search power vs. adapter alone
#   trim just the adapter portion, leaving the half-cut-site genomic bases intact
#   add tags to QNAME indicating channel and trim status
# input:
#   an unaligned SAM stream, specifically, a Dorado-compatible output stream, on STDIN
# output:
#   updated SAM stream on STDOUT
# notes:
#   this script processes unaligned ONT reads, i.e., never reverse complemented, so first bases correspond to the 5' adapter
#   all reads are expected to have 5' adapters, but some may be low quality bases, truncated, or chimeras split by Dorado
#   3' adapters are always shorter, often truncated, and not always present for incompletely sequenced fragments
#   most detected adapters are expected to occur at RE sites in ligFree (but not tagFree) libraries

# load dependencies
our $script = "trim_adapters";
our $error  = "$script error";
our ($matchScore, $mismatchPenalty, $gapOpenPenalty, $gapExtensionPenalty) = 
    (1,           -1,               -1.5,            -1); # penalties lowered from defaults of -1.5, -2.5, -1
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general smith_waterman);

# environment variables
fillEnvVar(\our $BLUNT_RE_TABLE,         'BLUNT_RE_TABLE'); # ligFree only supports blunt REs; 5' REs in tagFree aren't used for endpoint trimming
fillEnvVar(\our $ENZYME_NAME,            'ENZYME_NAME');
fillEnvVar(\our $CHECK_ENDPOINT_RE_MATCH,'CHECK_ENDPOINT_RE_MATCH');
fillEnvVar(\our $ADAPTER_SEQUENCE,       'ADAPTER_SEQUENCE');

# set operating parameters
my $MIN_SCORE = 10; # 3 RE bases + 1 A-tail + 6 adapter bases
my $isCustomEnzyme = (uc($ENZYME_NAME) eq "CUSTOM");
my $includeREHalfSite = (!$isCustomEnzyme and $CHECK_ENDPOINT_RE_MATCH);

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
    MIN_INSERT_SIZE => 250
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
    $cutSiteLeftHalf or die "unrecognized enzyme; must be a blunt cutter in shared/modules/REs/blunt_enzymes.csv, or 'custom'\n";
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
        my @read = split("\t", $read, 12);
        length($read[SEQ]) >= MIN_INSERT_SIZE or next; # reject reads too small to even search for adapaters
        my ($found5, $trimLen5) = find_adapter($adapter5, substr($read[SEQ], 0, SEARCH_SPACE_5), FORCE_QRY_END,   5);
        my ($found3, $trimLen3) = find_adapter($adapter3, substr($read[SEQ], -SEARCH_SPACE_3),   FORCE_QRY_START, 3);
        if($found3){ # order is important
            $read[SEQ]  = substr($read[SEQ],  0, -$trimLen3);
            $read[QUAL] = substr($read[QUAL], 0, -$trimLen3);
        }
        if($found5){ # order is important
            $read[SEQ]  = substr($read[SEQ],  $trimLen5);
            $read[QUAL] = substr($read[QUAL], $trimLen5);
        }
        $read[QNAME] = join(":", 
            $read[QNAME], # QNAME exits as QNAME:channel:trim5:trim3
            $read[TAGS] =~ m/ch:i:(\d+)/ ? $1 : 0,
            $found5 ? $trimLen5 : 0, 
            $found3 ? $trimLen3 : 0
        );
        print join("\t", @read);
    }
}

# execute an efficient adapter search
# searches force alignment to the end of the adapter closest to the adapter-gDNA transition for specificity (also the highest quality bases)
sub find_adapter {
    my ($adapter, $readSegment, $forcedEnd, $end) = @_;
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
            # allow this on the 3' end only - 5' ends are expected to match an RE half-site
            # } elsif($end == 3) {
        } elsif ($includeREHalfSite) { # this would be redundant for tagFree, Cas9, etc.
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

# 414239e5-00ca-44cd-8276-489a8cdf3c12    4       *       0       0       *       *       0       0       
# ATGTGCCCTCTACTGGTTCAGTTAGTATTGCTATCTTACAAACTCTTTGAGACCATGAATACCAGCACAGATTACTACATCCAGCAAAATTTTCAACTACTGTAGATGAAGAAAACAAAACATTTCATGACAAAGTCAAATTTGAACATTATCCAACAATATAGCTGTACAGAAAATACTAAAAGGAAAATTCCAACCCAAGGAAGTTAACTACATCCATGAAAACATAAGCAATAAATAATCCTGCACCAGTATCACCCAAAGAAGGGATGTTTATTAATATATTTCAATATCAATTGACTAAATTCTCTAGTAAATAGGCACAGGCTAATAGAATGTATATAATCACAGAATCAGTATGAATATATGAAACATAACTCAATATTAAAGATAGATATTACCTCAGAATAAAGAGTTAGGAACAGATTTTTCCAAGCAAATGGACACACAAAGCAAACTGATATAACTGTTTTAAAATCTAACAAAGTAGCTTTCATATAAAATTAATCAAAAGATGGAGAAGTTCCTACTAATTAAAGAAAAAAATCTACCAAGATAATGTTTCAATTCATACCATCTATGTCCAAAATGCAAGAGAACCCACATTTGTAAAAGAAACATTTGTAAAGCTTATATAACATATTGAAACTTACATATTAATAATGGGAGATTTCAATACTTGATTCTGACCAGACAAAAAGTCAAAGAGAAATAAGAGAACTAGCAGTCAC      
# #$$%$&*++**)*+((((',/++(&&&(**+<<<:<<<=ACC>>===BBC><<<<@>?@A@><<99::;???@?A?@?@DDB?=ACCCOLIIJFFEEDCAA???>>@@@E@?@@PJDEEDDFCA=;;==DA@BA;;;;<FFJHEDDDEGFED<;<<<CGGGEE@@@>?AAABCBCIHHA?<<;:3>::<HFGC>>>=<<<<;<<<;<=;<=<=>>?>>>?@ECHNLKJJHDAAABCLMLKIH@>><;<:9:77899?;::;<CBCA:8>88=DDB@BFJJGGEDDDBABBBCDDBBDDA@>>>@ACHCB>==<=?@ADD:9888;;87889AB;;;;<ABBBFFIJEEBB???A@D@???;===>CCDCBCCBDCBB??>>>BBHEEFG@>==;<<<<EA>>>=>>?@@@DE@?>==@=<<=<>@@?==:978777;;<;DB>===>A???A<<<<=DDBABBBFDBB@@@BCFEIIMDBA@?AAC>==99::;BGFEEFEIJIKGFHGBA?:5*'%%$$$$$%&&&),/266885558<<=GMJJIBA@?>AA@ABAADBBCBABBAA====<>>>=?AA>@>>=?@CGA???=><<<=;<<<<ABCFD///0>?SODDDEEA????A????AGFFFDDCCCDBBBBCFEDCCCCCEFIHIFFF@>=<<>ABBCA@?@@??@CECCBB@AB@@AAABDFGEFFCMKIHGGHGD>=::;200.......-- 
# qs:i:20 du:f:1.8524     ns:i:9262       ts:i:10 mx:i:4  ch:i:1092       
# st:Z:2024-04-25T20:38:44.242+00:00   rn:i:51659      
# fn:Z:PAO64638_pass_006c84b6_9f549c95_1140.pod5  sm:f:-810.172   
# sd:f:0.00798761 sv:Z:pa dx:i:0  RG:Z:9f549c95bcae9f37c4365b25d9a7ef0166b3a379_dna_r10.4.1_e8.2_400bps_sup@v4.3.0

use strict;
use warnings;

# action:
#   PENDING: if requested, reject sequences that aligned better to another species in a zoo
#   order alignments across each read from 5' to 3'
#   record metadata on number of ref/read bases in alignment
# input:
#   name-sorted SAM on STDIN
# output: 
#   name-sorted SITE_SAM on STDOUT, ordered from read 5' to 3' end

# initialize reporting
our $script = "order_alignments";
our $error  = "$script error";
my ($nSequences) = (0) x 10;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
# fillEnvVar(\our $ZOO_FILTER_LIST,   'ZOO_FILTER_LIST');

# constants
use constant {
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
    # -------------
    FALSE   => 0,
    TRUE    => 1,
};
my $nullMetadata = join("\t", 
    1, # BLOCK_N
    (0) x (TARGET_CLASS - BLOCK_N) # SITE_INDEX1_1..TARGET_CLASS
);

# load any available zoo rejections (sequence that aligned better to another species)
my %zooRejections;
# if(-f $ZOO_FILTER_LIST){
#     open my $zooH, "-|", "zcat $ZOO_FILTER_LIST" or die "could not open $ZOO_FILTER_LIST: $!\n";
#     while (my $qName = <$zooH>){
#         chomp $qName;
#         $zooRejections{$qName}++;
#     }
#     close $zooH;
# }
my $nZooRejections = scalar(keys %zooRejections);

# parse input SAM
my ($prevQName, @alns, @outAlns);
while(my $aln = <STDIN>){
    chomp $aln;
    my @aln = split("\t", $aln, SPLIT_TO_TAGS);

    # execute zoo rejection
    $zooRejections{$aln[QNAME]} and next;

    # handle sequences as a set of alignments
    if($prevQName and $prevQName ne $aln[QNAME]){
        processQName();
        @alns = @outAlns = ();
    }
    push @alns, \@aln;
    $prevQName = $aln[QNAME];
}
processQName(); 

# print summary information
# printCount(commify($nZooRejections),'nZooRejections',   'sequences rejected by zoo');
printCount(commify($nSequences),    'nSequences',       'analyzed sequences');

# process QNAME alignment groups
sub processQName {
    $nSequences++;
    if(@alns == 1){
        @outAlns = ($alns[0]);
    } else {
        @outAlns = sort { getQueryStart0(@{$a}[FLAG, CIGAR]) <=> getQueryStart0(@{$b}[FLAG, CIGAR]) } @alns;
    }
    print join("\n", map {
        my $ch = $outAlns[$_][TAGS] =~ m/ch:i:(\d+)/ ? $1 : 0;
        my $tl = $outAlns[$_][TAGS] =~ m/tl:Z:(\S+)/ ? $1 : "0,0";
        my $de = $outAlns[$_][TAGS] =~ m/de:f:(\S+)/ ? $1 : 0;
        my $hv = $outAlns[$_][TAGS] =~ m/hv:i:(\d+)/ ? $1 : 0;
        my $cs = $outAlns[$_][TAGS] =~ m/cs:Z:(\S+)/ ? $1 : "*";
        join("\t", 
            @{$outAlns[$_]}[QNAME..CIGAR], 
            $ch, $tl, $de, $hv,
            getRefSpan($outAlns[$_][CIGAR]),     # N_REF_BASES
            getAlignedSize($outAlns[$_][CIGAR]), # N_READ_BASES
            $nullMetadata, 
            @{$outAlns[$_]}[SEQ..QUAL],
            $cs
        );
    } 0..$#outAlns), "\n";
}

1;

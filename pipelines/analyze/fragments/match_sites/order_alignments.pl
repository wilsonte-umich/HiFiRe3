use strict;
use warnings;

# action:
#   PENDING: for alignments generated externally, add null extensions to QNAME
#   if requested, reject sequences that aligned better to another species in a zoo
#   for paired reads:
#       check read gaps for evidence of chimeras/SVs and if found convert to two orphans
#           if a gap is split, adjust xf:i: tag so that event bits reflect only each orphaned read
#       change the readN of read2 orphans created here or upstream to read1
#   add QNAME:splitGapReadN extension to all sequences reflecting paired gap chimera splitting, where
#       QNAME:0 means sequence was not split at a gap (unpaired platforms, quality/unmapped orphans, and proper pairs)
#       QNAME:N=1|2 means sequence was split at a gap, where N is the original readN for each orphan
#   order alignments across each read from 5' to 3'
#   record metadata on number of ref/read bases in alignment
#   note that read gaps are uninformative for SVs when using RE ends because we
#       cannot perform RE-based SV error correction if a junction was not sequenced
#       cannot use TLEN distribution to assess likelihood of non-sequenced insertions/deletions
# input:
#   name-sorted SAM on STDIN
# output: 
#   name-sorted SITE_SAM on STDOUT, ordered from read 5' to 3' end, where:
#       QNAME is appended with :splitGapReadN
#       SITE_SAM comprises fields (many with default values when leaving this script)
#           QNAME to CIGAR      same as SAM
#           DE_TAG              minimap2 de:f: tag = fractional gap-compressed per-base sequence divergence
#           CS_TAG              minimap2 cz:Z: tag = short-format but complete reference alignment encoding, i.e., difference string
#           XF_TAG              custom xf:i: tag = bit encoding of EVENT_HAS_VARIANT, EVENT_HAS_SV, etc
#           XH_TAG              custom xh:i: tag = bit encoding of matching haplotypes as determined by alignment
#           N_REF_BASES         number of reference bases spanned by the alignment
#           N_READ_BASES        number of read bases comprising the alignment
#           BLOCK_N             reference vs. read traversal block, numbered sequentially across each read prior to quality truncation
#           SITE_INDEX1_1       1-indexed site on RNAME for 5'-most alignment end pos1 in read order, signed by the direction from pos1 to SITE_POS1_1
#           SITE_POS1_1         1-indexed RE site coordinate on RNAME for SITE_INDEX1_1
#           SITE_HAPS_1         bit-encoded haplotypes that have a site matching SITE_INDEX1_1
#           SITE_DIST_1         signed distance in bp from strand-adjusted (but NOT clip-adjusted) pos1 to SITE_POS1_1
#           SITE_INDEX1_2       same as above for 3'-most alignment end pos1 in read order
#           SITE_POS1_2         "
#           SITE_HAPS_2         "
#           SITE_DIST_2         "
#           SEQ_SITE_INDEX1_2   unsigned SITE_INDEX1_2 for the actual or projected sequence 3' end
#           SEQ_SITE_POS1_2     1-indexed RE site coordinate on RNAME for SEQ_SITE_INDEX1_2
#           SEQ_SITE_HAPS_2     bit-encoded haplotypes that have a site matching SEQ_SITE_INDEX1_2
#           IS_END_TO_END       boolean 0|1 flag, true if the sequence 3' end reached a RE site to complete a fragment
#           READ_HAS_JXN        boolean 0|1 flag, true if the read retains an SV junction AFTER chimera splitting WITHOUT haplotype overide
#           TARGET_CLASS        bitwise encoding of whether the alignment and read overlapped a target region (1-bit = alignment, 2-bit = read)
#           SEQ to QUAL         same as SAM, often * on non-variant reads
#           (no TAGS, all metadata found in custom SITE_SAM columns)

# initialize reporting
our $script = "order_alignments";
our $error  = "$script error";
my ($nSequences, $nPairGaps, $nSplitGaps, $nRead2Orphans) = (0) x 10;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
# fillEnvVar(\our $IS_USER_BAM,       'IS_USER_BAM');
fillEnvVar(\our $READ_PAIR_TYPE,    'READ_PAIR_TYPE');
fillEnvVar(\our $FIXED_READ_LENGTH, 'FIXED_READ_LENGTH');
fillEnvVar(\our $ZOO_FILTER_LIST,   'ZOO_FILTER_LIST');
my $MAX_INSERT_SIZE = getMaxInsertSize();
$FIXED_READ_LENGTH = $FIXED_READ_LENGTH || 0;
my $MAX_GAP_LENGTH = $READ_PAIR_TYPE eq "paired" ? $MAX_INSERT_SIZE - 2 * $FIXED_READ_LENGTH : 0;
# my $nullQNameExtension = ":".join(":", 
#     0, # channel = external reads assume non-ONT platform
#     0, # trim5
#     0, # trim3
#     0, # isMerged = assume external reads were never merged
#     0, # nRead1
#     0, # nRead2
#     # splitGapReadN = added by this script
# );

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
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
    #-------------
    ALN_HAS_SNV         => 1, # variant flag bits, xf:i
    READ_HAS_SNV        => 2,
    READ_HAS_SV         => 4,
    EVENT_HAS_SNV       => 8,
    EVENT_HAS_SV        => 16,
    EVENT_HAS_VARIANT   => 32,
    # -------------
    FALSE   => 0,
    TRUE    => 1,
};
my $readFlagBits    = READ_HAS_SV + READ_HAS_SNV;
my $alnReadFlagBits = READ_HAS_SV + READ_HAS_SNV + ALN_HAS_SNV;
my $nullMetadata = join("\t", 
    1, # BLOCK_N
    (0) x (TARGET_CLASS - BLOCK_N) # SITE_INDEX1_1..TARGET_CLASS
);

# load any available zoo rejections (sequence that aligned better to another species)
my %zooRejections;
if(-f $ZOO_FILTER_LIST){
    open my $zooH, "-|", "zcat $ZOO_FILTER_LIST" or die "could not open $ZOO_FILTER_LIST: $!\n";
    while (my $qName = <$zooH>){
        chomp $qName;
        $zooRejections{$qName}++;
    }
    close $zooH;
}
my $nZooRejections = scalar(keys %zooRejections);

# parse input SAM
my ($prevQName, @alns, @outAlns, @outSplitNs);
my $nSeqAlns = 0;
while(my $aln = <STDIN>){
    chomp $aln;
    my @aln = split("\t", $aln, SPLIT_TO_TAGS);
    # $IS_USER_BAM and $aln[QNAME] .= $nullQNameExtension;

    # execute zoo rejection
    $zooRejections{$aln[QNAME]} and next;

    # handle sequences as a set of alignments
    if($prevQName and $prevQName ne $aln[QNAME]){
        processQName();
        @alns = @outAlns = @outSplitNs = ();
        $nSeqAlns = 0;
    }
    my $readN = ($aln[FLAG] & _IS_PAIRED and $aln[FLAG] & _SECOND_IN_PAIR) ? READ2 : READ1;
    push @{$alns[$readN]}, \@aln;
    $prevQName = $aln[QNAME];
    $nSeqAlns++;
}
processQName(); 

# print summary information
printCount(commify($nZooRejections),'nZooRejections',   'sequences rejected by zoo');
printCount(commify($nSequences),    'nSequences',       'input sequences after zoo rejection');
printCount(commify($nPairGaps),     'nPairGaps',        'paired reads with gaps');
printCount(commify($nSplitGaps),    'nSplitGaps',       'anomalous gaps split at inferred unverifiable junctions');
printCount(commify($nRead2Orphans), 'nRead2Orphans',    'sequence with read2-only orphans, e.g., read1 unmapped');

# process QNAME alignment groups
sub processQName {
    $nSequences++;

    # paired reads
    if($alns[READ1] and $alns[READ2]){
        $nPairGaps++;
        foreach my $readN(READ1, READ2){ # query order here is defined by each read individually
            @{$alns[$readN]} > 1 and @{$alns[$readN]} = sort { getQueryStart0(@{$a}[FLAG, CIGAR]) <=> getQueryStart0(@{$b}[FLAG, CIGAR]) } @{$alns[$readN]};
        }

        # anomalous paired gaps, split to two orphans
        if(isAnomalousGap(
            $alns[READ1][$#{$alns[READ1]}], # innermost alignments on each read define the gap
            $alns[READ2][$#{$alns[READ2]}]
        )){
            $nSplitGaps++;
            foreach my $aln(@{$alns[READ2]}){  # when pair was split as an anomalous gap ...
                $$aln[FLAG] ^= _FIRST_IN_PAIR; # ... reassign newly orphaned read2 as read1
                $$aln[FLAG] ^= _SECOND_IN_PAIR;
            }
            @outSplitNs = (
                (1) x @{$alns[READ1]},         # ... and retain the original readN as a flag that QNAME was split
                (2) x @{$alns[READ2]}
            );

        # proper paired gaps, retain as pair
        } else {
            @outSplitNs = (0) x $nSeqAlns;
        }
        @outAlns = (@{$alns[READ1]}, @{$alns[READ2]}); # each read individually handled in 5' to 3' order

    # single reads, merged read pairs, or read1 quality/unmapped orphans
    } elsif($alns[READ1]){
        if($nSeqAlns == 1){
            @outAlns = ($alns[READ1][0]);
        } else {
            @outAlns = sort { getQueryStart0(@{$a}[FLAG, CIGAR]) <=> getQueryStart0(@{$b}[FLAG, CIGAR]) } @{$alns[READ1]};
        }
        @outSplitNs = (0) x @outAlns;

    # read2 quality/unmapped orphans
    } else {
        $nRead2Orphans++;
        foreach my $aln(@{$alns[READ2]}){
            $$aln[FLAG] ^= _FIRST_IN_PAIR;  # when read1 was lost, reassign orphaned read2 as sole read1
            $$aln[FLAG] ^= _SECOND_IN_PAIR; # thus, orphans always lose read2 from a pair
        }
        if($nSeqAlns == 1){
            @outAlns = ($alns[READ2][0]);
        } else {
            @outAlns = sort { getQueryStart0(@{$a}[FLAG, CIGAR]) <=> getQueryStart0(@{$b}[FLAG, CIGAR]) } @{$alns[READ2]};
        }
        @outSplitNs = (0) x @outAlns;
    }

    # print modified alignments and interleaved alignment metadata
    print join("\n", map {
        $outAlns[$_][QNAME] .= ":$outSplitNs[$_]";
        my ($de, $cs, $xf, $xh) = ($outAlns[$_][TAGS] =~ m/de:f:(\S+)\tcs:Z:(\S+)\txf:i:(\d+)\txh:i:(\d)/);
        $outSplitNs[$_] and $xf = adjustXfTag($xf);
        join("\t", 
            @{$outAlns[$_]}[QNAME..CIGAR], 
            $de, $cs, $xf, $xh,
            getRefSpan($outAlns[$_][CIGAR]),     # N_REF_BASES
            getAlignedSize($outAlns[$_][CIGAR]), # N_READ_BASES
            $nullMetadata, 
            @{$outAlns[$_]}[SEQ..QUAL]
        );
    } 0..$#outAlns), "\n";
}

# return isAnomalous, i.e., 1 if there appears to be a junction in the gap, 0 otherwise
sub isAnomalousGap { 
    my ($aln1, $aln2) = @_;
    $$aln1[RNAME] ne $$aln2[RNAME] and return TRUE; # translocation
    my $strand1 = ($$aln1[FLAG] & _REVERSE);
    my $strand2 = ($$aln2[FLAG] & _REVERSE);
    $strand1 == $strand2 and return TRUE;           # inversion
    my $refPos1_1 = $strand1 ? $$aln1[POS1] : getEnd(@{$aln1}[POS1, CIGAR]); # i.e., the inner node positions flanking the gap
    my $refPos1_2 = $strand2 ? $$aln2[POS1] : getEnd(@{$aln2}[POS1, CIGAR]); 
    my $gapLength = $strand1 ? 
        $refPos1_1 - $refPos1_2 - 1 :
        $refPos1_2 - $refPos1_1 - 1;
    (
        $gapLength < -$FIXED_READ_LENGTH or      # duplication
        $gapLength >  $MAX_GAP_LENGTH            # deletion
    ) and return TRUE;
    FALSE; 
}

# when paired reads were split to two orphans, adjust xf:i: tag to reflect each read individually
sub adjustXfTag {
    my ($xf) = @_;
    $xf or return $xf;
    my $readVariants = ($xf & $readFlagBits);
    ($xf & $alnReadFlagBits) + 
    ($readVariants << 2) + 
    ($readVariants ? EVENT_HAS_VARIANT : 0);
}

1;

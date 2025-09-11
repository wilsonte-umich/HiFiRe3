use strict;
use warnings;

# SV parsing functions shared by RE and non-RE fragment parsers

# working variables
use vars qw($SAMPLE_PATHS_PREFIX_WRK $paddedChromIndex1 $MAX_JXN_BASES %chromIndex);

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
    #-------------
    NULL_ENTRY      => "*",
    NULL_JUNCTION   => "*,*",
    NULL_JXN_BASES  => "#", # the value used when junction bases were not processed, but alignment offset is present
    #-------------
    NODE1 => 0,
    NODE2 => 1,
    #-------------
    NODE_CHROM_INDEX1   => 0,
    NODE_REF_POS1       => 1,
    NODE_STRAND_INDEX0  => 2,
    #-------------
    JXN_ALN_OFFSET   => 0,
    JXN_JXN_BASES    => 1,
    JXN_JXN_TYPE     => 2, # added after indexing
    #-------------
    TYPE_PROPER         => 0, # junction/edge types
    TYPE_DELETION       => 1,
    TYPE_DUPLICATION    => 2,
    TYPE_INVERSION      => 3,
    TYPE_TRANSLOCATION  => 4,
    #-------------
    UNKNOWN_HAPLOTYPE   => 0,
    HAPLOTYPE1          => 1,
    HAPLOTYPE2          => 2,
    REFERENCE_BIT       => 1,
    HAPLOTYPE1_BIT      => 2,
    HAPLOTYPE2_BIT      => 4,
    #-------------
    channel             => 0, # incoming qName extensions
    trim5               => 1,
    trim3               => 2,
    isMerged            => 3, # true=2 for legacy reasons
    nRead1              => 4,
    nRead2              => 5,
    splitGapReadN       => 6,
    readN               => 7,
    N_QNAME_EXTENSIONS  => 8,
};

# support conversion of haplotype bit-encoding to the final best haplotype
our $genotypeBits = HAPLOTYPE1_BIT + HAPLOTYPE2_BIT;
our @matchingHaps = (UNKNOWN_HAPLOTYPE);                             # 0
$matchingHaps[REFERENCE_BIT]                  = UNKNOWN_HAPLOTYPE;  # 1
$matchingHaps[HAPLOTYPE1_BIT]                 = HAPLOTYPE1;         # 2
$matchingHaps[HAPLOTYPE1_BIT + REFERENCE_BIT] = HAPLOTYPE1;         # 3
$matchingHaps[HAPLOTYPE2_BIT]                 = HAPLOTYPE2;         # 4
$matchingHaps[HAPLOTYPE2_BIT + REFERENCE_BIT] = HAPLOTYPE2;         # 5
$matchingHaps[$genotypeBits]                  = UNKNOWN_HAPLOTYPE;  # 6
$matchingHaps[$genotypeBits  + REFERENCE_BIT] = UNKNOWN_HAPLOTYPE;  # 7

# open/close the required chromosome-level output files
our ($jxnSrcH, $uniqJxnH, $firstAlnH, $distalAlnH, $alnSizeH, $readPathH);
sub openFileHandles {
    open $jxnSrcH,    "|-", "gzip -c > $SAMPLE_PATHS_PREFIX_WRK.junction_sources.$paddedChromIndex1.txt.gz"  or die "could not open jxnSrcH: $!\n";
    open $uniqJxnH,   "|-", "gzip -c > $SAMPLE_PATHS_PREFIX_WRK.unique_junctions.$paddedChromIndex1.txt.gz"  or die "could not open uniqJxnH: $!\n";
    open $firstAlnH,  "|-", "gzip -c > $SAMPLE_PATHS_PREFIX_WRK.first_alignments.$paddedChromIndex1.txt.gz"  or die "could not open firstAlnH: $!\n";
    open $distalAlnH, "|-", "gzip -c > $SAMPLE_PATHS_PREFIX_WRK.distal_alignments.$paddedChromIndex1.txt.gz" or die "could not open distalAlnH: $!\n";
    open $alnSizeH,   "|-", "gzip -c > $SAMPLE_PATHS_PREFIX_WRK.alignment_sizes.$paddedChromIndex1.txt.gz"   or die "could not open alnSizeH: $!\n";
    open $readPathH,  "|-", "gzip -c > $SAMPLE_PATHS_PREFIX_WRK.read_paths.$paddedChromIndex1.txt.gz"        or die "could not open readPathH: $!\n";
}
sub closeFileHandles {
    close $jxnSrcH;
    close $uniqJxnH;
    close $firstAlnH;
    close $distalAlnH;
    close $readPathH;
}

# functions to describe paths
sub getJunction {
    my ($aln1, $aln2) = @_; # alignments to read 5' and 3' of a junction, respectively

    # cannot describe junction when SEQ was dropped by `align genotype`, i.e., haplotype-adjusted READ_HAS_SV is FALSE
    $$aln1[S_SEQ] eq NULL_ENTRY and return NULL_JUNCTION;

    # otherwise, determine the offset/overlap between the two alignments
    # alignmentOffset is positive for insertions, negative for microhomologies, 0 for blunt
    my $clip1 = 
        ($$aln1[S_FLAG] & _REVERSE) ? 
            getLeftClip( $$aln1[S_CIGAR]) :
            getRightClip($$aln1[S_CIGAR]);
    my $clip2 = 
        ($$aln2[S_FLAG] & _REVERSE) ? 
            getRightClip($$aln2[S_CIGAR]) :
            getLeftClip( $$aln2[S_CIGAR]);
    my $alignmentOffset = ($clip1 + $clip2) - length($$aln1[S_SEQ]); 

    # then determine the read bases corresponding to the insertion or microhomology
    # note that inserted bases must be provided by the read, not reference alignments
    # JXN_BASES leave in the orientation of the original read, thus,
    # will have two distinct rc'ed values when then same junction was sequenced in opposite orientations
    my $jxnBases = NULL_ENTRY; # *, a blunt joint
    # if(abs($alignmentOffset) > $MAX_JXN_BASES){
    #     $jxnBases = NULL_JXN_BASES; # large number of bases from platforms with low base accuracy are rarely useful and ~always have errors
    # } els
    if($alignmentOffset > 0){
        $jxnBases = ($$aln1[S_FLAG] & _REVERSE) ? 
            getRc( substr(substr($$aln1[S_SEQ], 0, $clip1), -$alignmentOffset) ) :
                   substr(substr($$aln1[S_SEQ], -$clip1), 0, $alignmentOffset);
    } elsif($alignmentOffset < 0) {
        $jxnBases = ($$aln1[S_FLAG] & _REVERSE) ? 
            getRc( substr(substr($$aln1[S_SEQ], 0, ($clip1 - $alignmentOffset)),   $alignmentOffset) ):
                   substr(substr($$aln1[S_SEQ], -($clip1 - $alignmentOffset)), 0, -$alignmentOffset);
    }
    join(",", $alignmentOffset, $jxnBases);
}

# enforce strict ordering of junction nodes for downstream junction aggregation
sub orientJunction {
    my ($node1, $node2, $junction, $junctionTypes, $isSplit) = @_;
    my @nodes = 
        $isSplit ? 
        (
            $node1,
            $node2
        ) : 
        (
            [ split(",", $node1) ],
            [ split(",", $node2) ]
        );
    my @junction = split(",", $junction);
    my $junctionType = getJunctionType(@nodes); # nodes before re-ordering
    $junctionTypes and push @$junctionTypes, $junctionType;

    # enforce strict ordering of junction pairs
    # translocations always sort low chromI to high chromI, i.e., by refPos1 along the composite genome
    # nearly all same-chrom junctions sort low refPos1 to high refPos1
    # only very rare inversions with identical refPos1 sort by strand, largely negligible
    # as a consequence, the canonical ordering of junctions is:
    #   deletions     sort to +/+, i.e., 1--->~~~2---> as canonical, <---2~~~<---1 gets flipped to it
    #                                        *   *                       *   *
    #   duplications  sort to -/-, i.e., <---1~~~<---2 as canonical, 2--->~~~1---> gets flipped to it
    #                                    *           *               *           *
    #   FF inversions sort to +/-, i.e., 1--->~~~<---2 as canonical, <---2~~~1---> gets flipped to it
    #                                        *       *                   *       *
    #   RR inversions sort to -/+, i.e., <---1~~~2---> as canonical, 2--->~~~<---1 gets flipped to it
    #                                    *       *                   *       *
    # the above patterns apply the same to translocations, where ~~|~~ simply crosses (a) chromosome boundary(ies)
    #     i.e., translocations follow either a del, dup, inv-FF or inv-RR pattern recognizable by oriented strand0 values
    # FF inversions can be viewed as "left anchored",  where read1 matches the reference genome orientation on the left side
    # RR inversions can be viewed as "right anchored", where read2 matches the reference genome orientation on the right side
    my @nodeI = sort { 
        $nodes[$a][NODE_CHROM_INDEX1]  <=> $nodes[$b][NODE_CHROM_INDEX1] or 
        $nodes[$a][NODE_REF_POS1]      <=> $nodes[$b][NODE_REF_POS1] or
        $nodes[$a][NODE_STRAND_INDEX0] <=> $nodes[$b][NODE_STRAND_INDEX0]
    } NODE1, NODE2;

    # when junction is re-oriented relative to source read ...
    if($nodeI[NODE1] == NODE2){ 
        # ... flip node strands
        $nodes[NODE1][NODE_STRAND_INDEX0] = (1 - $nodes[NODE1][NODE_STRAND_INDEX0]);
        $nodes[NODE2][NODE_STRAND_INDEX0] = (1 - $nodes[NODE2][NODE_STRAND_INDEX0]);
        # ... reverse complement junction bases
        # continuing above, now expect junctions sequenced in different orientations to have the same JXN_BASES in canonical read orientation
        rc(\$junction[JXN_JXN_BASES]);
        # note that S_SEQ/JXN_SEQ, i.e., the sequences of different reads matching the junction, are not rc'ed here
        # thus reads of multiply sequenced junctions may be in two different orientations as revealed by wasInverted
    }

    # return the ordered junction
    (
        @{$nodes[$nodeI[NODE1]]}, # ordered node1 = chrom pos strand
        @{$nodes[$nodeI[NODE2]]}, # ordered node2 = chrom pos strand
        $junctionType, 
        @junction,    # alnOffset jxnBases
        $nodeI[NODE1] # strandIndex0, i.e., wasInverted
    )
}
sub getJunctionType {
    my ($node1, $node2) = @_; # nodes before re-ordering
    $$node1[NODE_CHROM_INDEX1]  != $$node2[NODE_CHROM_INDEX1]  and return TYPE_TRANSLOCATION;
    $$node1[NODE_STRAND_INDEX0] != $$node2[NODE_STRAND_INDEX0] and return TYPE_INVERSION;
    if($$node1[NODE_STRAND_INDEX0] == TOP_STRAND){
        $$node1[NODE_REF_POS1]  <  $$node2[NODE_REF_POS1]      and return TYPE_DELETION;
    } else {
        $$node2[NODE_REF_POS1]  <  $$node1[NODE_REF_POS1]      and return TYPE_DELETION;
    }
    return TYPE_DUPLICATION;
}

# functions to describe paths
sub getRefPos1_5 {
    my ($aln) = @_;
    ($$aln[S_FLAG] & _REVERSE) ? 
        getEnd(@$aln[S_POS1, S_CIGAR]) :
        $$aln[S_POS1];
}
sub getNode1 {
    my ($aln) = @_;
    my $strandIndex0 = ($$aln[S_FLAG] & _REVERSE) ? BOTTOM_STRAND : TOP_STRAND;
    join(",", 
        $chromIndex{$$aln[S_RNAME]},
        $strandIndex0 == TOP_STRAND ? 
            $$aln[S_POS1] :
            getEnd(@$aln[S_POS1, S_CIGAR]),
        $strandIndex0
    )
}
sub getNode2 {
    my ($aln) = @_;
    my $strandIndex0 = ($$aln[S_FLAG] & _REVERSE) ? BOTTOM_STRAND : TOP_STRAND;
    join(",", 
        $chromIndex{$$aln[S_RNAME]},
        $strandIndex0 == TOP_STRAND ? 
            getEnd(@$aln[S_POS1, S_CIGAR]):
            $$aln[S_POS1],
        $strandIndex0
    )
}

# parse pseudo-reference sequence contig for one SV read across all alignments
sub getPseudoRefSeq {
    my ($alns) = @_;
    my $maxAlnsI = $#{$alns};
    join("", map { 
        my $i = $_;
        my $seq = getRefSeq2(
            $$alns[$i][S_RNAME], 
            $$alns[$i][S_POS1], 
            getEnd($$alns[$i][S_POS1], $$alns[$i][S_CIGAR]), 
            $$alns[$i][S_FLAG] & _REVERSE
        );
        if($i < $maxAlnsI){
            # must get junction anew, since we require analysis of even junctions filtered by BLOCK_N above
            my ($alignmentOffset, $jxnBases) = split(",", getJunction($$alns[$i], $$alns[$i + 1]));
            if($alignmentOffset < 0){ # microhomology bases contributed by aln2 for each junction, trimmed from aln1
                $seq = substr($seq, 0, length($seq) + $alignmentOffset);
            } elsif($alignmentOffset > 0) {
                $seq .= $jxnBases; # these and only these junction insertion bases in pseudo reference contig come from the read itself
            }
        }
        $seq
    } 0..$maxAlnsI);
}

# outer endpoint positions and ONT channels for duplicate purging
sub getOuterPos5 {
    my ($aln) = @_;
    my $strandIndex0 = ($$aln[S_FLAG] & _REVERSE) ? BOTTOM_STRAND : TOP_STRAND;
    $strandIndex0 == TOP_STRAND ? 
        $$aln[S_POS1] :
        getEnd(@$aln[S_POS1, S_CIGAR])
}
sub getOuterPos3 {
    my ($aln) = @_;
    my $strandIndex0 = ($$aln[S_FLAG] & _REVERSE) ? BOTTOM_STRAND : TOP_STRAND;
    $strandIndex0 == TOP_STRAND ? 
        getEnd(@$aln[S_POS1, S_CIGAR]) :
        $$aln[S_POS1]
}
sub getONTChannel {
    my ($aln) = @_;
    my @qName = split(":", $$aln[S_QNAME]);
    my @extensions = splice(@qName, -N_QNAME_EXTENSIONS);
    $extensions[channel]
}
sub getOuterEndpoints {
    my ($chromIndex1_5, $outerPos1_5, $chromIndex1_3, $outerPos1_3, $channel) = @_;
    my $oe1 = join("\t", $chromIndex1_5, $outerPos1_5);
    my $oe2 = join("\t", $chromIndex1_3, $outerPos1_3);
    ($oe1 cmp $oe2) <= 0 ? # yields the same key ordering for both strands
        (join("\t", $oe1, $oe2, $channel), 0) : # report whether we reordered endpoints
        (join("\t", $oe2, $oe1, $channel), 1);
}
sub getOuterEndpoints_lig_free {
    my ($chromIndex1_5, $chromIndex1_3) = @_;
    $chromIndex1_5 <= $chromIndex1_3 ? # yields the same key ordering for both strands
        (join("\t", $chromIndex1_5, $chromIndex1_3), 0) : # report whether we reordered endpoints
        (join("\t", $chromIndex1_3, $chromIndex1_5), 1);
}

1;

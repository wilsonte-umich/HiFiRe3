use strict;
use warnings;

# actions:
#   analyze singleton junctions to find SVs that match another SV path or a local non-SV path sufficiently to reject them
#   each singleton junction sequence is re-aligned against:
#       - its own pseudo-reference contig, i.e., the null hypothesis that minimap2 got it right
#       - the reference sequence at breakpoint 1, extended as appropriate for alternative continuous local alignment at aln1
#       - the reference sequence at breakpoint 2, extended as appropriate for alternative continuous local alignment at aln2
#       - the pseudo-reference contig of other junctions matching at one breakpoint, to ask if the singleton matches another validated junction
#   in all cases, the query junction sequence is aligned to (pseudo-)reference sequences (not to other reads)
#   for consistency, all target sequences are:
#       - ~the same reference length, mathching the length of the query junction sequence
#       - padded similarly to allow account for the impact of indels in achieving complete alignments
#   which means that each query singleton junction re-alignment must be handled individually
#   criteria for accepting an alternative alignment, thus masking the singleton junction:
#       - alternative alignment score >= 90% of the null hypothesis alignment score
#       - for local alternative alignment:
#           - the side being added must account for >= 50% of the length of that same side of the null hypothesis alignment
# input:
#   junction stream on STDIN
# output:
#   junction stream on STDOUT with new column UJXN_HAS_ALT_ALIGNMENT set to TRUE|FALSE

# initialize reporting
our $script = "cross_reference";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general faidx);

# environment variables
fillEnvVar(\our $N_CPU, 'N_CPU');
fillEnvVar(\our $GENOME_FASTA_SHM,    'GENOME_FASTA_SHM');
fillEnvVar(\our $TMP_FILE_PREFIX_SHM, 'TMP_FILE_PREFIX_SHM');
fillEnvVar(\our $ALIGNMENT_MODE,      'ALIGNMENT_MODE');
fillEnvVar(\our $BANDWIDTH,           'BANDWIDTH');
fillEnvVar(\our $GROUP_BREAKPOINT_DISTANCE, 'GROUP_BREAKPOINT_DISTANCE');
fillEnvVar(\our $CROSS_REFERENCE_PADDING,   'CROSS_REFERENCE_PADDING');

# constants
use constant {
    OJXN_CHROM_INDEX1_1     => 0, # input columns
    OJXN_REF_POS1_1         => 1,
    OJXN_STRAND_INDEX0_1    => 2,
    OJXN_CHROM_INDEX1_2     => 3,
    OJXN_REF_POS1_2         => 4,
    OJXN_STRAND_INDEX0_2    => 5,
    OJXN_JXN_TYPE           => 6,
    UJXN_N_OBSERVED         => 7,
    UJXN_ALN_OFFSET         => 8,
    UJXN_JXN_BASES          => 9,
    UJXN_PATHS              => 10,
    UJXN_N_PATH_JUNCTIONS   => 11,
    UJXN_READ_HAS_SV        => 12,
    UJXN_JSRC_MAPQ          => 13, 
    UJXN_JSRC_DE_TAG        => 14,
    UJXN_JSRC_SITE_DIST     => 15,
    UJXN_JSRC_QNAMES        => 16, 
    UJXN_JSRC_SEQS          => 17,
    UJXN_JSRC_QUALS         => 18,
    UJXN_JSRC_CIGARS        => 19,
    UJXN_JSRC_ORIENTATIONS  => 20,
    UJXN_HAS_ALT_ALIGNMENT  => 21, # added by this script
    #-----------------
    BKPT_JXN_ID        => 0,
    BKPT_CHROM_INDEX1  => 1,
    BKPT_REF_POS1      => 2,
    BKPT_DIRECTION     => 3,
    #-----------------
    QUERY_JXN_SEQ   => "qry",
    NULL_HYPOTHESIS => "nll",
    ALT_LOCAL_1     => "loc-1",
    ALT_LOCAL_2     => "loc-2",
    ALT_OTHER_JXN   => "jxn",
    #-----------------
    TOP_STRAND    => 0,
    BOTTOM_STRAND => 1,
    #-----------------
    ALN_1 => 0,
    ALN_2 => 1,
    #-----------------
    SIDE_1 => 1,
    SIDE_2 => 2,
    #-----------------
    LEFTWARD  => 0,
    RIGHTWARD => 1,
    #-----------------
    FALSE => 0,
    TRUE  => 1,
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
    SPLIT_TO_TAGS => 12,
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
    RESULT_N_BASES   => 0,
    RESULT_ALN_SCORE => 1,
    #-------------
    ALN_SCORE_THRESHOLD   => 0.9,
    ALN_LOC_LEN_THRESHOLD => 0.5,
};

# initialize the genome
print STDERR "collecting chrom index\n";
use vars qw(%revChromIndex);
setCanonicalChroms();
loadFaidx($GENOME_FASTA_SHM);
open our $faH, "<", $GENOME_FASTA_SHM or die "could not open $GENOME_FASTA_SHM: $!\n";

# initalize minimap2 realignment
my $fqFile = "$TMP_FILE_PREFIX_SHM.query.fq";
my $prFile = "$TMP_FILE_PREFIX_SHM.pseudo_ref.fa";
my $minimap2Command = "minimap2 -ax $ALIGNMENT_MODE -t $N_CPU -Y --secondary=yes -N 10 -r $BANDWIDTH $prFile $fqFile 2>/dev/null";

# load junctions, collecting all breakpoints into a single list
my (
    $jxnId, $prevLocusKey, $prevRefPos1, @jxnIds, @jxnTargets, @jxns, @breakpoints, 
    $nFlankBases, $nQryBases_1, $nQryBases_2, $expNBases, $qryWasInverted
) = (0);
while (my $jxn = <STDIN>){
    $jxnId++;
    chomp $jxn;
    my @jxn = split("\t", $jxn);
    $jxn[UJXN_HAS_ALT_ALIGNMENT] = FALSE; # initialize as no alternative alignment, overridden as needed beloe
    $jxns[$jxnId] = \@jxn;
    # convert strands to directions so the chrom/pos/direction is the same regardless of bkpt pairing and reordering
    push @breakpoints, [
        $jxnId, 
        @jxn[OJXN_CHROM_INDEX1_1, OJXN_REF_POS1_1], 
        OJXN_STRAND_INDEX0_1 == TOP_STRAND ? LEFTWARD : RIGHTWARD
    ];
    push @breakpoints, [
        $jxnId, 
        @jxn[OJXN_CHROM_INDEX1_2, OJXN_REF_POS1_2], 
        OJXN_STRAND_INDEX0_2 == TOP_STRAND ? RIGHTWARD : LEFTWARD
    ];
}

# collect the set of junctions with a fuzzy breakpoint match on _either_ side of each junction
# this is not the same as prior group_nodes, which matched on _both_ ordered breakpoints
# note that bkpt1 is compared to both bkpt1 and bkpt2 of other junctions
# since a given alignment could end up as bkpt1 or bkpt2 depending on breakpoint reordering
foreach my $bkpt(sort {
    $$a[BKPT_CHROM_INDEX1] <=> $$b[BKPT_CHROM_INDEX1] or
    $$a[BKPT_DIRECTION]    <=> $$b[BKPT_DIRECTION] or
    $$a[BKPT_REF_POS1]     <=> $$b[BKPT_REF_POS1]
} @breakpoints){
    my $locusKey = join(":", @{$bkpt}[BKPT_CHROM_INDEX1, BKPT_DIRECTION]);
    if($prevLocusKey and 
    (
        $prevLocusKey ne $locusKey or 
        $prevRefPos1 < $$bkpt[BKPT_REF_POS1] - $GROUP_BREAKPOINT_DISTANCE # establish position breaks
    )
    ){
        processBreakpointGroup();
        @jxnIds = ();
    }
    push @jxnIds, $$bkpt[BKPT_JXN_ID];
    $prevLocusKey = $locusKey;
    $prevRefPos1 = $$bkpt[BKPT_REF_POS1];
}
processBreakpointGroup();
@breakpoints = ();
sub processBreakpointGroup {
    foreach my $qryJxnId(@jxnIds){
        foreach my $tgtJxnId(@jxnIds){
            $jxnTargets[$qryJxnId]{$tgtJxnId}++; # includes de-duplicated self-to-self pairs 
        }                                        # and both reciprocal pairs, e.g., 2-3 and 3-2 as different test cases
    }
}

# run the realignment of each junction to the alternative pseudo-reference contigs
foreach my $qryJxnId(1..$#jxnTargets){ # there is no junction ID 0

    # open minimap2 temporary file handles
    open my $fqH, ">", $fqFile or die "$error: could not open $fqFile : $!\n";
    open my $prH, ">", $prFile or die "$error: could not open $prFile : $!\n"; # $faH variable name used for the reference genome above

    # only evaluate singleton junctions as query junctions
    # junctions with multiple observations are considered sufficiently validated
    my $qryJxn = $jxns[$qryJxnId];
    if($$qryJxn[UJXN_N_OBSERVED] > 1){
        print join("\t", @$qryJxn), "\n";
        next;
    }

    # extract just the read bases that comprise the singleton query junction alignment
    # despite the singleton restriction, UJXN_JSRC_SEQS, _QUALS, _CIGARS and _ORIENTATIONS can have multiple entries due to de-duplication
    my $qrySeq  = getFirstReadElement($qryJxn, UJXN_JSRC_SEQS); # oriented as per the original 5'-most read alignment
    my $qryQual = getFirstReadElement($qryJxn, UJXN_JSRC_QUALS);
    my @cigars = split("::", getFirstReadElement($qryJxn, UJXN_JSRC_CIGARS)); # in 5'-3' read order
    $qryWasInverted = getFirstReadElement($qryJxn, UJXN_JSRC_ORIENTATIONS);
    my ($strandIndex0_1, $strandIndex0_2) = $qryWasInverted ?
        (1 - $$qryJxn[OJXN_STRAND_INDEX0_2], 1 - $$qryJxn[OJXN_STRAND_INDEX0_1]) :
        (    $$qryJxn[OJXN_STRAND_INDEX0_1],     $$qryJxn[OJXN_STRAND_INDEX0_2]);
    if($strandIndex0_1 == BOTTOM_STRAND){ # UJXN_JSRC_SEQ is oriented as per aln1; restore all to original read 5'-3' order
        rc(\$qrySeq);
        $qryQual = reverse($qryQual);
    }
    my $qryStart0 = getQueryStart0( # thus, the first query base of the first, 5'-most junction alignment
        $strandIndex0_1 == BOTTOM_STRAND ? _REVERSE : 0, 
        $cigars[ALN_1]
    );
    my $qryEnd1 = getQueryEnd1( # thus, the last query base of the second, 3'-most junction alignment
        $strandIndex0_2 == BOTTOM_STRAND ? _REVERSE : 0, 
        $cigars[ALN_2], 
        length($qrySeq)
    );
    my $qryJxnSeq  = substr($qrySeq,  $qryStart0, $qryEnd1 - $qryStart0); # implicitly includes UJXN_JXN_BASES
    my $qryJxnQual = substr($qryQual, $qryStart0, $qryEnd1 - $qryStart0);
    my $fqName = join("_", QUERY_JXN_SEQ, $qryJxnId);
    print $fqH '@'."$fqName\n$qryJxnSeq\n+\n$qryJxnQual\n";

    # collect all competing pseudo-reference contigs for re-alignment of the single query junction
    # for simplicity, use the same number of bases on each side of the junction for each contig
    $nQryBases_1 = getAlignedSize($cigars[ALN_1]);
    $nQryBases_2 = getAlignedSize($cigars[ALN_2]);
    $nFlankBases = max($nQryBases_1, $nQryBases_2) + max(0, $$qryJxn[UJXN_ALN_OFFSET]) + $CROSS_REFERENCE_PADDING;
    foreach my $tgtJxnId(keys %{$jxnTargets[$qryJxnId]}){

        # add the pseudo-reference contig for the junction to its original split alignment
        # this is the null hypothesis that the initial minimap2 alignment was correct
        if($qryJxnId == $tgtJxnId){
            my $prName = NULL_HYPOTHESIS;
            my $prSeq = getJxnPseudoReference($qryJxn);
            print $prH ">$prName\n$prSeq\n";

            # add the alternative pseudo-reference contig assuming that kept aln1 should be expanded as a single locus alignment
            $prName = ALT_LOCAL_1;
            $prSeq = getLocalPseudoReference($qryJxn, SIDE_1);
            print $prH ">$prName\n$prSeq\n";

            # add the alternative pseudo-reference contig assuming that kept aln2 should be expanded as a single locus alignment
            $prName = ALT_LOCAL_2;
            $prSeq = getLocalPseudoReference($qryJxn, SIDE_2);
            print $prH ">$prName\n$prSeq\n";

        # add the pseudo-reference contig for the junction to another junction with a single breakpoint match
        } else {
            my $tgtJxn = $jxns[$tgtJxnId];
            my $prName = join("_", ALT_OTHER_JXN, $tgtJxnId);
            my $prSeq = getJxnPseudoReference($tgtJxn);
            print $prH ">$prName\n$prSeq\n";
        }
    }
    close $fqH;
    close $prH;

    # perform the realignment of the query junction to the alternative pseudo-reference contigs
    open my $mm2H, "-|", "$minimap2Command" or die "$error: could not open minimap2 stream : $!\n";

    # parse realignments, collecting null hypothesis and alternative alignment metadata
    my ($prevQname, %results);
    $expNBases = $nQryBases_1 + $nQryBases_2 + $$qryJxn[UJXN_ALN_OFFSET];
    while (my $aln = <$mm2H>){
        $aln =~ m/^@/ and next; # skip header lines
        chomp $aln;
        my @aln = split("\t", $aln, SPLIT_TO_TAGS);
        ($aln[FLAG] & _SUPPLEMENTAL) and next;
        $aln[TAGS] =~ m/AS:i:(\S+)/;
        my $alnScore = $1;
        my $prType = (split("_", $aln[RNAME]))[0];
        if($results{$prType}){ # the second jxn type encountered
            $results{$prType}[RESULT_ALN_SCORE] >= $alnScore and next; # already better than this alignment
        } 
        $results{$prType} = [getAlignedSize($aln[CIGAR]), $alnScore];
        $prevQname = $aln[QNAME];
    }

    # determine if a junction has an alternative nearly as good as (or better than) the null hypothesis
    if($results{nll}){
        if(
            getJxnResults(ALT_LOCAL_1,   $results{'loc-1'}, $results{'nll'}) or
            getJxnResults(ALT_LOCAL_2,   $results{'loc-2'}, $results{'nll'}) or
            getJxnResults(ALT_OTHER_JXN, $results{'jxn'},   $results{'nll'})
        ){
            $$qryJxn[UJXN_HAS_ALT_ALIGNMENT] = TRUE;
        }
    } else { # unexpected code path, but if there was no null hypothesis alignment, there must have been a good alternative alignment
        $$qryJxn[UJXN_HAS_ALT_ALIGNMENT] = TRUE;
    }
    print join("\t", @$qryJxn), "\n";
    close $mm2H;
}
sub getFirstReadElement {
    my ($jxn, $col) = @_;
    (split(",", $$jxn[$col]))[0];
}


# finish up
close $faH;

# get a portion of the reference genome for assembling a pseudo-reference contig
sub getPseudoRefSegment {
    my ($chromIndex1, $refPos1, $strandIndex0, $side) = @_;
    my $direction = $side == SIDE_1 ? 
        ($strandIndex0 == TOP_STRAND ? LEFTWARD : RIGHTWARD) : 
        ($strandIndex0 == TOP_STRAND ? RIGHTWARD : LEFTWARD);
    getRefSeq2( # this function re-orders positions as needed
        $revChromIndex{$chromIndex1}, 
        $refPos1, 
        $direction == RIGHTWARD ? 
            $refPos1 + $nFlankBases - 1 : 
            $refPos1 - $nFlankBases + 1, 
        $strandIndex0
    );
}

# assemble a pseudo-reference contig for a specific target junction, sized according to the query junction
# the target and query junctions can be the same or different
sub getJxnPseudoReference {
    my ($tgtJxn) = @_;
    my $seg1 = getPseudoRefSegment(@{$tgtJxn}[OJXN_CHROM_INDEX1_1, OJXN_REF_POS1_1, OJXN_STRAND_INDEX0_1], SIDE_1);
    my $seg2 = getPseudoRefSegment(@{$tgtJxn}[OJXN_CHROM_INDEX1_2, OJXN_REF_POS1_2, OJXN_STRAND_INDEX0_2], SIDE_2);
    my $tgtJxnSeq;
    if($$tgtJxn[UJXN_ALN_OFFSET] < 0){
        $seg1 = substr($seg1, 0, length($seg1) + $$tgtJxn[UJXN_ALN_OFFSET]); # remove microhomology bases from one side
        $tgtJxnSeq = $seg1.$seg2;
    } elsif($$tgtJxn[UJXN_ALN_OFFSET] == 0){
        $tgtJxnSeq = $seg1.$seg2;
    } else {
        $tgtJxnSeq = $seg1.$$tgtJxn[UJXN_JXN_BASES].$seg2; # UJXN_JXN_BASES already oriented relative to reoriented target aln1
    }
    $tgtJxnSeq;
}

# assemble a pseudo-reference contig for extended local alignment on one side of the query junction
sub getLocalPseudoReference {
    my ($qryJxn, $keptSide) = @_;
    my ($seg1, $seg2);
    if($keptSide == SIDE_1){
        # the kept alignment reference
        $seg1 = getPseudoRefSegment(@{$qryJxn}[OJXN_CHROM_INDEX1_1, OJXN_REF_POS1_1, OJXN_STRAND_INDEX0_1], SIDE_1);
        # the locally extended alternative alignment reference
        $seg2 = getPseudoRefSegment(@{$qryJxn}[OJXN_CHROM_INDEX1_1, OJXN_REF_POS1_1, OJXN_STRAND_INDEX0_1], SIDE_2);
    } else {
        # the locally extended alternative alignment reference
        $seg1 = getPseudoRefSegment(@{$qryJxn}[OJXN_CHROM_INDEX1_2, OJXN_REF_POS1_2, OJXN_STRAND_INDEX0_2], SIDE_1);
        # the kept alignment reference
        $seg2 = getPseudoRefSegment(@{$qryJxn}[OJXN_CHROM_INDEX1_2, OJXN_REF_POS1_2, OJXN_STRAND_INDEX0_2], SIDE_2);
    }
    $seg1.$seg2;
}

# assess whether an alternative alignment is nearly as good as (or better) than the null hypothesis
# return TRUE is this is a valid alternative alignment
sub getJxnResults {
    my ($prType, $result, $nllResult) = @_;
    $result or return FALSE;

    # first check the minimap2 alignment score
    my $passed = ($$result[RESULT_ALN_SCORE] >= $$nllResult[RESULT_ALN_SCORE] * ALN_SCORE_THRESHOLD);

    # for local alignments, further check that the alternative alignment accounts for
    # a sufficient fraction of the bases in the local extension, i.e., the side opposite the kept alignment
    if($passed and $prType eq ALT_LOCAL_1){
        my $keptSideLen = $qryWasInverted ? $nQryBases_2 : $nQryBases_1;
        my $altSideLen  = $qryWasInverted ? $nQryBases_1 : $nQryBases_2;
        $passed = ($$result[RESULT_N_BASES] - $keptSideLen >=  $altSideLen * ALN_LOC_LEN_THRESHOLD);
    }
    if($passed and $prType eq ALT_LOCAL_2){
        my $keptSideLen = $qryWasInverted ? $nQryBases_1 : $nQryBases_2;
        my $altSideLen  = $qryWasInverted ? $nQryBases_2 : $nQryBases_1;
        $passed = ($$result[RESULT_N_BASES] - $keptSideLen >=  $altSideLen * ALN_LOC_LEN_THRESHOLD);
    }
    $passed ? TRUE : FALSE;
}

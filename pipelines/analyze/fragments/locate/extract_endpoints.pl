use strict;
use warnings;

# action:
#     extract all sequence endpoints expected to match RE sites, i.e.:
#       5' endpoints of all sequences, i.e., reads
#       3' endpoints from 3'-adapter-trimmed single-reads or merged read pairs, i.e., molecules sequenced end-to-end
#           *>>>>>>>>>>>o------  (truncated/untrimmed single-read, 3' end not used)
#           *>>>>>>>>>>>>>>>>>>* (merged read pair, 3'-adapter-trimmed ONT or other end-to-end single read)
#           *>>>>>>------<<<<<<* (unmerged paired reads)
#     where "endpoint" means the outermost bases of the sequenced event, ignoring internal supplementary alignments
#           *aln1//aln2//aln3* (aln2 is ignored by this script)
# expects:
#     source $MODULES_DIR/genome/set_genome_vars.sh
#     source $MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/library/set_library_vars.sh
# output: 
#     headerless table of observed endpoints with columns:
#           chrom,sitePos1,nObserved
#     where:
#       chrom is the chromosome name per setCanonicalChroms()
#       sitePos1 is the 1-referenced genome position adjusted to infer at least one endpoint's candidate blunt RE site
#                  |*======   where | = RE cleaved bond, + = refPos1, * = sitePos1, side = R
#           ======+|*         and side = L
#                  |*--+===   where clipped (-) but not trimmed bases imply that RE site is offset from refPos1
#           ====+--|*         small clips more likely reflect base differences in contiguous alignments, large clips may represent SV junctions
#       nObserved is the count of all unique endpoints that nominated chrom,sitePos1
#     ordered by chrom, but not by sitePos1 (since will be reordered anyway by tabulate_endpoints.R)
#     all sitePos1 are for blunt REs, since read endpoint analysis only applies to ligation libraries

# initialize reporting
our $action = "extract_endpoints";
my ($nAlns, $nSequences, $nReads, $nEndpoints, $nHighQual, $nPositions) = (0) x 10;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\my $SEQUENCING_PLATFORM,    'SEQUENCING_PLATFORM');
fillEnvVar(\my $IS_END_TO_END_READ,     'IS_END_TO_END_READ');
fillEnvVar(\my $MIN_MAPQ,               'MIN_MAPQ');
fillEnvVar(\my $CLIP_TOLERANCE,         'CLIP_TOLERANCE');
fillEnvVar(\my $NAME_BAM_FILE,          'NAME_BAM_FILE');

# set platform-specific parameters
my $isEndToEndPlatform = ($IS_END_TO_END_READ eq "TRUE");
my $isONT = ($SEQUENCING_PLATFORM eq "ONT");
my $ontAdapterLen5 = 34; # empirically determined median adapter trim lengths 

# assistance for adjusting refPos1 to sitePos1 based on clips
#        *           *           *     * = sitePos1
#   ----|xx5-----3xx/xx5-----3xx|----  as nodes are numbered for alignments
#   ----|xx3-----5xx/xx3-----5xx|----
my %clipDir = (
    "5:+" => -1,
    "5:-" =>  1,
    "3:+" =>  1,
    "3:-" => -1
);

# initialize the genome
use vars qw(%chromIndex %revChromIndex @canonicalChroms);
setCanonicalChroms();

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
    QSTART0 => 12,
    #-------------
    SPLIT_TO_TAGS => 12,
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
    TRIM5 => 0, 
    TRIM3 => 1,
    #-------------
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
};

# run all reads in all input PAF files
print STDERR "extracting sequence endpoints\n";
my ($prevQName, @alns, $readN, @counts) = (0);
$| = 1;
open my $inH, "-|", "samtools view $NAME_BAM_FILE" or die "$action error: could not open $NAME_BAM_FILE: $!\n";
while(my $aln = <$inH>){
    $nAlns++;
    my @aln = split("\t", $aln, SPLIT_TO_TAGS);
    $aln[TAGS] = "\t".$aln[TAGS]; # facilitate regex searches
    if($prevQName and $prevQName ne $aln[QNAME]){
        $nSequences++;
        foreach my $readN (READ1..READ2) { 
            $alns[$readN] and processRead($readN) 
        }
        @alns = ();
    }
    $aln[QSTART0] = getQueryStart0(@aln[FLAG, CIGAR]);
    $readN = ($aln[FLAG] & _IS_PAIRED and $aln[FLAG] & _SECOND_IN_PAIR) ? READ2 : READ1;
    push @{$alns[$readN]}, \@aln;
    $prevQName = $aln[QNAME];
}
close $inH;

# process the last sequence
$nSequences++;
foreach my $readN (READ1..READ2) { 
    $alns[$readN] and processRead($readN) 
}

# report all unique observed positions with counts
foreach my $chrom(@canonicalChroms){
    my $rNameIndex = $chromIndex{$chrom};
    my $nObserved = $counts[$rNameIndex] or next;
    foreach my $sitePos1(keys %$nObserved){
        print join("\t", $chrom, $sitePos1, $$nObserved{$sitePos1}), "\n";
        $nPositions++;
    }
}

# report tallies
printCount(commify($nAlns),      'nAlns',      'alignments');
printCount(commify($nSequences), 'nSequences', 'sequences, i.e., reads');
printCount(commify($nReads),     'nReads',     'reads');
printCount(commify($nEndpoints), 'nEndpoints', 'informative endpoints');
printCount(commify($nHighQual),  'nHighQual',  "candidate endpoints with MAPQ >= $MIN_MAPQ and clip <= $CLIP_TOLERANCE");
printCount(commify($nPositions), 'nPositions', "genome positions nominated as RE sites");

# process a read over one or two outermost alignments
sub processRead {
    my ($readN) = @_;
    $nReads++;

    # sort alignments by query start position
    my @alns_ = sort { $$a[QSTART0] <=> $$b[QSTART0] } @{$alns[$readN]};

    # get ONT adapter trims
    my @trims = (0, 0);
    $isONT and $alns_[0][TAGS] =~ /\ttl:Z:(\d+),(\d+)/ and @trims = ($1, $2);

    # always attempt to process the 5' end of all reads
    processEndpoint($alns_[0], \@trims, 5);

    # attempt to process read 3' ends if informative, i.e., inferred to be from end-to-end read
    # NOT necessarily on the same chromosome or part of the same contig as the 5' end
    $readN == READ2 and return;
    if(
        $IS_END_TO_END_READ or        # PacBio or other platforms that guarantee end-to-end reads
        ($isONT and $trims[TRIM3]) or # complete ONT reads, must be trimmed at 3' end to be confident end-to-end
        $alns_[0][TAGS] =~ /\tfm:Z:2/ # merged paired reads
    ){
        processEndpoint($alns_[$#alns_], \@trims, 3);
    }
}

# check an informative endpoint for sufficient alignment quality
sub processEndpoint {
    my ($aln, $trims, $end) = @_;
    $nEndpoints++;
    ($chromIndex{$$aln[RNAME]} and $$aln[MAPQ] >= $MIN_MAPQ) or return;
    my $strand = ($$aln[FLAG] & _REVERSE) ? "-" : "+";
    my $clipLen = 
        $end == 5 ? 
        $$aln[QSTART0] : (
            $strand eq "+" ? 
            getRightClip($$aln[CIGAR]) :
            getLeftClip($$aln[CIGAR])
        );
    $clipLen <= $CLIP_TOLERANCE or return;
    $nHighQual++;

    # parse the combination of endpoint rPos1, strand, and clip to predicted sitePos1
    #        *           *           *     * = sitePos1
    #   ----|xx5-----3xx/xx5-----3xx|----  as nodes are numbered for alignments
    #   ----|xx3-----5xx/xx3-----5xx|----  
    my $nodeKey = join(":", $end, $strand);
    my $sitePos1 = 
        $end == 5 ? 
        ($strand eq "+" ? $$aln[POS1] : getEnd($$aln[POS1], $$aln[CIGAR]) + 1) + $clipDir{$nodeKey} * $clipLen :
        ($strand eq "+" ? getEnd($$aln[POS1], $$aln[CIGAR]) + 1 : $$aln[POS1]) + $clipDir{$nodeKey} * $clipLen;

    # use a typical trim value to adjust clips on ONT trim failures, to account for adapter bases still present
    # only applicable to 5' ends here, untrimmed 3' ends were filtered above as uninformative
    $isONT and $end == 5 and $$trims[TRIM5] == 0 and $sitePos1 -= $clipDir{$nodeKey} * $ontAdapterLen5;

    # keep a tally of all genome positions nominated as RE sites
    $counts[$chromIndex{$$aln[RNAME]}]{$sitePos1}++;
}

use strict;
use warnings;

# scope:
#   this script is applied to non-RE-cleaved fragments, i.e., a non-agFree library,
#   to skip non-applicable actions in discard_unmatched_sequences.pl
# action:
#   set projection fields and end-to-end flag
# input:
#   name-sorted SITE_SAM on STDIN
# output: 
#   name-sorted SITE_SAM on STDOUT
#   where this script:
#       updates SEQ_SITE_POS1_2 and IS_END_TO_END when known
#       appends readN to QNAME

# initialize reporting
our $script = "set_end_to_end";
our $error  = "$script error";
my ($nSequences) = (0) x 10;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\our $SEQUENCING_PLATFORM,    'SEQUENCING_PLATFORM');
fillEnvVar(\our $READ_PAIR_TYPE,         'READ_PAIR_TYPE');
fillEnvVar(\our $IS_END_TO_END_READ,     'IS_END_TO_END_READ');

# set platform-specific parameters
my $isONT = $SEQUENCING_PLATFORM eq "ONT";
my $isPairedReadPlatform = ($READ_PAIR_TYPE eq "paired");
my $isEndToEndPlatform = ($IS_END_TO_END_READ eq "TRUE");

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
    # -------------
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
    # -------------
    FALSE  => 0,
    TRUE   => 1,
};

# process alignments one at a time
my (@alns, $prevQName);
while(my $aln = <STDIN>){
    my @aln = split("\t", $aln, SPLIT_TO_IS_END_TO_END);
    if($prevQName and $prevQName ne $aln[S_QNAME]){
        processQName();
        @alns = ();
    }
    my $readN = ($aln[S_FLAG] & _IS_PAIRED and $aln[S_FLAG] & _SECOND_IN_PAIR) ? READ2 : READ1;
    push @{$alns[$readN]}, \@aln;
    $prevQName = $aln[S_QNAME];
}
processQName();

# print summary information
printCount(commify($nSequences), 'nSequences', 'input sequences');

# process sequences by QNAME
sub processQName {
    $nSequences++;

    # determine whether the reads together sequenced both end of the original insert
    my $aln5 = $alns[READ1][0];
    my @qName = split(":", $$aln5[S_QNAME]);
    my @extensions = splice(@qName, -N_QNAME_EXTENSIONS);
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
                my $aln3 = $alns[READ2] ? $alns[READ2][0] : $alns[READ1][$#{$alns[READ1]}];
                @$aln[SEQ_SITE_POS1_2] = 
                    $readN == READ1 ? 
                    (
                        $alns[READ2] ? 
                        $$aln3[SITE_POS1_1] :
                        $$aln3[SITE_POS1_2]
                    ) : 
                    @$aln5[SITE_POS1_1];
            }
            $$aln[IS_END_TO_END] = $isEndToEnd;
            print join("\t", @$aln);
        }
    }
}

use strict;
use warnings;

# action:
#   flag alignments, reads, and events for the presence or absence of SVs and SNV/indels
#   remove QUAL on an alignment when it has no SNVs
#   remove SEQ on an alignment when it has no SVs or SNVs (note: must retain SEQ if QUAL is retained)
#       except on ONT where SEQ and QUAL are always retained since ~all read have SNVs and QUAL is used for SV filtering
# input:
#   SAM from reference minimap2 alignment on STDIN
# output:
#   if requested for genotype, TMP_VARIANT_FASTQ, a file with reads from events with variants for realignment
#   modified SAM on STDOUT with custom tags added, at end, in order:
#       xf:i:   bit-encoded flag of event, read, and alignment variant status
#       xh:i:   bit-encoded flag of matching haplotypes (always leaves this script as 1, i.e., matched to reference)

# initialize reporting
our $script = "flag_variant_events";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);


# environment variables
fillEnvVar(\our $HAS_BASE_ACCURACY, 'HAS_BASE_ACCURACY');
fillEnvVar(\our $ALN_CPU, 'ALN_CPU');
my ($TMP_VARIANT_FASTQ) = @ARGV;

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
    HAS_SNV => 12, # temporary field not in output 
    #-------------
    SPLIT_TO_TAGS  => 12,
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
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
    #-------------
    ALN_HAS_SNV         => 1, # variant flag bits
    READ_HAS_SNV        => 2,
    READ_HAS_SV         => 4,
    EVENT_HAS_SNV       => 8,
    EVENT_HAS_SV        => 16,
    EVENT_HAS_VARIANT   => 32,
};
my $alnScoreDelim = "~~";  # don't use '/' or '|' as this delimiter

# open the output file
my $varH;
if($TMP_VARIANT_FASTQ){
    open $varH, "|-", "pigz -c --processes $ALN_CPU > $TMP_VARIANT_FASTQ" or die "could not open $TMP_VARIANT_FASTQ: $!\n";
}

# parse input SAM
my @hasSNV = (0, 0, 0);
my (@alns, @hasSV, $prevQName);
while(my $aln = <STDIN>){

    # pass header lines as is
    if($aln =~ m/^\@/){
        print $aln;
        next;
    }

    # process alignments
    chomp $aln;
    my @aln = split("\t", $aln, SPLIT_TO_TAGS);
    if($prevQName and $prevQName ne $aln[QNAME]){
        processQName();
        @hasSNV = (0, 0, 0);
        @alns = ();
    }
    my $readN = ($aln[FLAG] & _IS_PAIRED and $aln[FLAG] & _SECOND_IN_PAIR) ? READ2 : READ1;
    if($HAS_BASE_ACCURACY and                          # i.e., not working on ONT reads
       $aln[TAGS] and $aln[TAGS] =~ m/cs:Z:(\S+)/ and  # i.e., read is mapped
       $1 !~ m/^:\d+$/                                 # alignment is not perfect (ignoring terminal clips)
    ){
        $hasSNV[EVENT]  = EVENT_HAS_SNV; # or indel...
        $hasSNV[$readN] = READ_HAS_SNV;
        $aln[HAS_SNV]   = ALN_HAS_SNV;
    } else {
        $aln[HAS_SNV]   = 0; # ONT reads, unmapped reads, perfect alignments (ignoring terminal clips)
    }
    push @{$alns[$readN]}, \@aln;
    $prevQName = $aln[QNAME];
}
processQName(); 
$varH and close $varH;

# process QNAME alignment groups
sub processQName {

    # assess SVs; SVs in paired read gaps are NOT recorded
    @hasSV = (0, 0, 0);
                     @{$alns[READ1]} > 1 and $hasSV[READ1] = READ_HAS_SV; 
    $alns[READ2] and @{$alns[READ2]} > 1 and $hasSV[READ2] = READ_HAS_SV;
    ($hasSV[READ1] or $hasSV[READ2]) and $hasSV[EVENT] = EVENT_HAS_SV;
    my $eventHasVariant = ($hasSNV[EVENT] or $hasSV[EVENT]) ? EVENT_HAS_VARIANT : 0;
    my $eventFlag = $eventHasVariant + $hasSV[EVENT]  + $hasSNV[EVENT];

    # print variant events to FASTQ
    if($varH and $eventHasVariant){
        my $compositeReadAlnScores = join(
            $alnScoreDelim, 
            0, # event placeholder, unused
                           getCompositeReadAlnScore($alns[READ1]), 
            $alns[READ2] ? getCompositeReadAlnScore($alns[READ2]) : 0
        );
        foreach my $readN(READ1, READ2){
            $alns[$readN] or next;
            my $aln = $alns[$readN][0];
            print $varH join("\n",
                '@'.$$aln[QNAME]."$alnScoreDelim$compositeReadAlnScores", # append composite alignment scores for evaluating relative quality of haplotype alignments
                ($$aln[FLAG] & _REVERSE) ? getRc($$aln[SEQ])    : $$aln[SEQ],
                '+',
                ($$aln[FLAG] & _REVERSE) ? reverse($$aln[QUAL]) : $$aln[QUAL]
            ), "\n";
        }
    }

    # print all SAM to STDOUT
    # unmapped reads will never have alignment or read-level variant flags set
    # but event-level variant flags may be set on unmapped reads if a paired read had a variant
    foreach my $readN(READ1, READ2){
        $alns[$readN] or next;
        foreach my $aln(@{$alns[$readN]}){

            # for most platforms that can call SNVs:
            #   QUAL remains on alignments with SNVs
            #   SEQ  remains on alignments with SNVs or reads with SVs (SEQ must exist if QUAL exists)
            # UPDATE: QUAL is now always retained for downstream SV filtering and junction analysis
            if($HAS_BASE_ACCURACY){
                # $$aln[HAS_SNV] or $$aln[QUAL] = "*";
                unless($$aln[HAS_SNV] or $hasSV[$readN]){
                    $$aln[QUAL] = "*";
                    $$aln[SEQ]  = "*";
                }

            # for ONT (or other low-accuracy platforms):
            #   SEQ and QUAL remain on reads with SVs (SNVs NA)
            #   (QUAL used downstream when evaluating SV flanking alignments)
            } elsif(!$hasSV[$readN]) {
                $$aln[QUAL] = "*";
                $$aln[SEQ]  = "*";
            }
            
            print join("\t",
                @$aln[QNAME..TAGS],
                "xf:i:".($eventFlag + $hasSV[$readN] + $hasSNV[$readN] + $$aln[HAS_SNV]),
                "xh:i:1"
            ), "\n";
        }
    }
}

1;

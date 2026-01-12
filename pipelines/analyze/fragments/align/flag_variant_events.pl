use strict;
use warnings;

# action:
#   flag alignments, reads, and events for the presence or absence of SVs and SNV/indels
# input:
#   SAM from reference minimap2 alignment on STDIN
# output:
#   modified SAM on STDOUT with custom tags added at end:
#       hv:i:   bit-encoded flag of alignment and read variant status

# initialize reporting
our $script = "flag_variant_events";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);

# environment variables
fillEnvVar(\our $HAS_BASE_ACCURACY, 'HAS_BASE_ACCURACY');

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
    if($HAS_BASE_ACCURACY and                               # i.e., not working on ONT reads
       $aln[TAGS] and "\t$aln[TAGS]" =~ m/\tcs:Z:(\S+)/ and # i.e., read is mapped
       $1 !~ m/^:\d+$/                                      # alignment is not perfect (ignoring terminal clips)
    ){
        $aln[HAS_SNV]   = ALN_HAS_SNV; # or indel...
        $hasSNV[$readN] = READ_HAS_SNV;
        $hasSNV[EVENT]  = EVENT_HAS_SNV;
    } else {
        $aln[HAS_SNV]  = 0;
    }
    push @{$alns[$readN]}, \@aln;
    $prevQName = $aln[QNAME];
}
processQName(); 

# process QNAME alignment groups
sub processQName {

    # assess SVs; SVs in paired read gaps are NOT recorded
    @hasSV = (0, 0, 0);
                     @{$alns[READ1]} > 1 and $hasSV[READ1] = READ_HAS_SV; 
    $alns[READ2] and @{$alns[READ2]} > 1 and $hasSV[READ2] = READ_HAS_SV;
    ($hasSV[READ1] or $hasSV[READ2]) and $hasSV[EVENT] = EVENT_HAS_SV;
    my $eventHasVariant = ($hasSNV[EVENT] or $hasSV[EVENT]) ? EVENT_HAS_VARIANT : 0;
    my $eventFlag = $hasSNV[EVENT] + $hasSV[EVENT] + $eventHasVariant;

    # print all SAM to STDOUT
    foreach my $readN(READ1, READ2){
        $alns[$readN] or next;
        foreach my $aln(@{$alns[$readN]}){
            print join("\t",
                @$aln[QNAME..TAGS],
                "hv:i:".($$aln[HAS_SNV] + $hasSNV[$readN] + $hasSV[$readN] + $eventFlag)
            ), "\n";
        }
    }
}

1;

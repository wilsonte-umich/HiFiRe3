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
    ALN_HAS_SNV       => 1, # variant flag bits
    READ_HAS_SNV      => 2,
    READ_HAS_SV       => 4,
    READ_HAS_VARIANT  => 8,
};

# parse input SAM
my $readHasSNV = 0;
my (@alns, $prevQName);
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
        $readHasSNV = 0;
        @alns = ();
    }
    if($HAS_BASE_ACCURACY and                          # i.e., not working on ONT reads
       $aln[TAGS] and $aln[TAGS] =~ m/cs:Z:(\S+)/ and  # i.e., read is mapped
       $1 !~ m/^:\d+$/                                 # alignment is not perfect (ignoring terminal clips)
    ){
        $aln[HAS_SNV] = ALN_HAS_SNV;
        $readHasSNV   = READ_HAS_SNV;
    } else {
        $aln[HAS_SNV]  = 0;
    }
    push @alns, \@aln;
    $prevQName = $aln[QNAME];
}
processQName(); 

# process QNAME alignment groups
sub processQName {

    # assess SVs based on number of alignments per read
    my $readHasSV = @alns > 1 ? READ_HAS_SV : 0;
    my $readHasVariant = ($readHasSNV or $readHasSV) ? READ_HAS_VARIANT : 0;
    my $readFlag = $readHasSNV + $readHasSV + $readHasVariant;

    # print all SAM to STDOUT
    foreach my $aln(@alns){
        print join("\t",
            @$aln[QNAME..TAGS],
            "hv:i:".($$aln[HAS_SNV] + $readFlag)
        ), "\n";
    }
}

1;

use strict;
use warnings;

# action:
#   find the closest RE site to each alignment node/endpoint
# input:
#   name-sorted SITE_SAM on STDIN
# output: 
#   name-sorted SITE_SAM on STDOUT
#   where this script updates values SITE_INDEX1_1 to SITE_DIST_2

# initialize reporting
our $script = "match_nodes_to_sites";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);

# environment variables
fillEnvVar(\our $ENZYME_NAME,             'ENZYME_NAME');
fillEnvVar(\our $BLUNT_RE_TABLE,          'BLUNT_RE_TABLE');
fillEnvVar(\our $SITE_CHROM_DATA,         'SITE_CHROM_DATA');          # in order of usage: first access the chrom's data
fillEnvVar(\our $CLOSEST_SITE_LOOKUP_WRK, 'CLOSEST_SITE_LOOKUP_WRK');  # then find the closest site on the chrom to query pos1
fillEnvVar(\our $SITE_DATA_LOOKUP_WRK,    'SITE_DATA_LOOKUP_WRK');     # then acquire the position and matching haplotypes

# constants
use constant {
    chrom_              => 0, # site index table columns
    chromIndex_         => 1,
    nSites_             => 2,
    chromSize_          => 3,
    closestSiteOffset_  => 4, # one lookup offset per chrom
    siteDataOffset_     => 5,
    #-------------
    SEEK_SET               => 0,   # see perl seek documentation
    CLOSEST_SITE_PACKING   => "l", # signed integer (siteIndex1)
    BYTES_PER_CLOSEST_SITE => 4,
    SITE_DATA_PACKING      => "LSS", # unsigned integer (sitePos1), unsigned shorts (two flanking frag sizes)
    BYTES_PER_SITE_DATA    => 8,
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
    SPLIT_TO_SITE_DIST_2 => 20,
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
};
my @nullSite = (0, 0, 0);

# load the RE site metadata to properly handle blunt sites
# enzyme  strand  cut_site regex   offset  CpG_priority
# EcoRV   0       GATATC   GATATC  3       4  
# thus:
# blunt cutter: correction5 = 0
#            *      sitePos1
#        --3 5--    top strand alignment
#  EcoRV GAT^ATC
#        CTA^TAG
#        --5 3--    bottom strand alignment
open my $reH, "<", $BLUNT_RE_TABLE or die "could not open: $!\n";
my $header = <$reH>; # enzyme,strand,cut_site,regex,offset,CpG_priority,high_fidelity,site_length
while (my $line = <$reH>){
    my ($enzyme, $strand, $cut_site, $regex, $offset, $priority, $hiFidelity, $siteLength) = split(",", $line);
    $enzyme or next;
    $enzyme =~ s/\s//g;
    $enzyme eq $ENZYME_NAME or next;
    last;
}
close $reH;


# load the RE site lookup index
my (%chromData);
open my $inH, "<", $SITE_CHROM_DATA  or die "could not open file $SITE_CHROM_DATA : $!";
$header = <$inH>; # chrom,chromIndex,nSites,chromSize,closestSiteOffset,siteDataOffset
while (my $line = <$inH>){
    chomp $line;
    my @chrom = split("\t", $line);
    $chromData{$chrom[chrom_]} = \@chrom;
}
close $inH;

# open handles to the RE filtering binary site lookups
open my $siteIndexH, '<:', $CLOSEST_SITE_LOOKUP_WRK or die "could not open file: $!";
open my $siteDataH,  '<:', $SITE_DATA_LOOKUP_WRK    or die "could not open file: $!";
my ($siteIndexRaw, $siteDataRaw) = ("", "");

# process alignments one at a time
while(my $aln = <STDIN>){

    # get closest RE sites
    # adjust alignment positions for strand/end offset and RE overhang but NOT clip yet as clip handling differs for outer and inner nodes
    # node1 as recorded here is always before node2 in read order
    # processing map, where 
    #   | is where RE cleaves, defining outer endpoint nodes
    #   / is a junction,       defining inner junction nodes
    #        *       *       *     proper sitePos1 values
    #   ----|5-----3/5-----3|----  as nodes are numbered for alignments
    #   ----|3-----5/3-----5|----
    my @aln = split("\t", $aln, SPLIT_TO_SITE_DIST_2);
    if($aln[S_FLAG] & _REVERSE){
        @aln[SITE_INDEX1_1..SITE_DIST_2] = (
            findClosestSite($aln[S_RNAME], getEnd(@aln[S_POS1, S_CIGAR]) + 1),
            findClosestSite($aln[S_RNAME], $aln[S_POS1])
        );
    } else {
        @aln[SITE_INDEX1_1..SITE_DIST_2] = ( 
            findClosestSite($aln[S_RNAME], $aln[S_POS1]),
            findClosestSite($aln[S_RNAME], getEnd(@aln[S_POS1, S_CIGAR]) + 1)
        );
    }
    print join("\t", @aln);
}

# find and learn about the RE site closest to a clipped endpoint node as reported by the aligner
# this site:
#   may not be the most relevant sitePos1 yet for incomplete 3' ends (see later scripts for 3' projection)
#   is not yet adjusted for clips (see later scripts that reject or split sequences based on RE site matching)
sub findClosestSite {
    my ($chrom, $adjPos1) = @_;
    $chromData{$chrom} or return @nullSite;

    # get the signed index of the closest site on the chromosome, on any haplotype
    my $siteIndexOffset = ${$chromData{$chrom}}[closestSiteOffset_] + ($adjPos1 - 1) * BYTES_PER_CLOSEST_SITE;
    seek($siteIndexH, $siteIndexOffset, SEEK_SET) or die "could not seek in siteIndex lookup: $!\n";
    if(!read($siteIndexH, $siteIndexRaw, BYTES_PER_CLOSEST_SITE)){
        # print STDERR "$error: could not read from siteIndex lookup: $chrom, $adjPos1: $!\n";
        return @nullSite;
    }
    my $siteIndex1 = unpack(CLOSEST_SITE_PACKING, $siteIndexRaw);

    # get the corresponding sitePos1
    my $siteDataOffset = ${$chromData{$chrom}}[siteDataOffset_] + (abs($siteIndex1) - 1) * BYTES_PER_SITE_DATA;
    seek($siteDataH, $siteDataOffset, SEEK_SET) or die "could not seek in siteData lookup: $!\n";
    if(!read($siteDataH, $siteDataRaw, BYTES_PER_SITE_DATA)){
        # print STDERR "$error: could not read from siteIndex lookup: $chrom, $adjPos1: $!\n";
        return @nullSite;
    }
    my ($sitePos1) = unpack(SITE_DATA_PACKING, $siteDataRaw); # discard fragment sizes here

    # calculate the signed distance from adjPos1 to sitePos1 and return all values
    return ($siteIndex1, $sitePos1, $adjPos1 - $sitePos1);
}

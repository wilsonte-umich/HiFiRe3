use strict;
use warnings;

# actions:
#   find the closest RE site to each alignment node/endpoint
#   project 3' clipped ends to the next closest RE site

# environment variables
use vars qw(
    $SITE_CHROM_DATA_FILE
    $CLOSEST_SITE_LOOKUP_WRK
    $SITE_DATA_LOOKUP_WRK
); # $ENZYME_NAME $BLUNT_RE_TABLE

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
    S_POS1              => 3,
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    SITE_INDEX1_1       => 6,
    SITE_POS1_1         => 7,
    SITE_DIST_1         => 8,
    SITE_INDEX1_2       => 9,
    SITE_POS1_2         => 10,
    # ... unused fields ...
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
    P_SITE_INDEX1  => 0,
    P_SITE_POS1    => 1,
};
my @nullSite = (0, 0, 0);

# # load the RE site metadata to properly handle blunt sites
# # enzyme  strand  cut_site regex   offset  CpG_priority
# # EcoRV   0       GATATC   GATATC  3       4  
# # thus:
# # blunt cutter: correction5 = 0
# #            *      sitePos1
# #        --3 5--    top strand alignment
# #  EcoRV GAT^ATC
# #        CTA^TAG
# #        --5 3--    bottom strand alignment
# open my $reH, "<", $BLUNT_RE_TABLE or die "could not open: $!\n";
# my $header = <$reH>; # enzyme,strand,cut_site,regex,offset,CpG_priority,high_fidelity,site_length
# while (my $line = <$reH>){
#     my ($enzyme, $strand, $cut_site, $regex, $offset, $priority, $hiFidelity, $siteLength) = split(",", $line);
#     $enzyme or next;
#     $enzyme =~ s/\s//g;
#     $enzyme eq $ENZYME_NAME or next;
#     last;
# }
# close $reH;

# load the RE site lookup index
our (%chromData);
open my $inH, "<", $SITE_CHROM_DATA_FILE  or die "could not open file $SITE_CHROM_DATA_FILE : $!";
my $header = <$inH>; # chrom,chromIndex,nSites,chromSize,closestSiteOffset,siteDataOffset
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

# find and learn about the RE site closest to a clipped endpoint node as reported by the aligner
# this site:
#   may not be the most relevant sitePos1 yet for incomplete 3' ends (see later scripts for 3' projection)
#   is not yet adjusted for clips (see later steps that reject or split sequences based on RE site matching)
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

# examine 5' most and 3'-most read alignments for a prior site projection of their 3' ends
# if missing, project them now to the next RE sites
# record whether:
#     these projections are in a read with a junction, i.e., the nature of the 3' nodes being projected
#     each alignment and read overlapped a target region
# alignment projection is indicated to fill out RE fragments for plotting and insert size tallies when:
#     a sequence was a truncated, i.e., was an incomplete read on a platform with variable-length single reads
#     a sequence had an SV that disconnected the 5'-most and 3'-most alignments of the parent sequence
# end projections are used downstream during fragment indexing and path contruction
# where all outermost endpoints of all sequences are described as the best matching RE site and orientation
sub getProjection {
    my ($aln) = @_;
    my @proj = @{$aln}[SITE_INDEX1_2..SITE_POS1_2];
    if($$aln[S_FLAG] & _REVERSE){     # no sign on either, pretend an exact match
        $proj[P_SITE_INDEX1] < 0 and  # actual and projected site indices could be the same
        abs($proj[P_SITE_INDEX1]) > 1 and
        $proj[P_SITE_INDEX1] = abs($proj[P_SITE_INDEX1]) - 1;
    } else {
        $proj[P_SITE_INDEX1] > 0 and 
        $proj[P_SITE_INDEX1] < ${$chromData{$$aln[S_RNAME]}}[nSites_] and 
        $proj[P_SITE_INDEX1]++; 
    }
    $proj[P_SITE_INDEX1] < 0 and $proj[P_SITE_INDEX1] = abs($proj[P_SITE_INDEX1]);
    $proj[P_SITE_INDEX1] == abs($$aln[SITE_INDEX1_2]) or
        $proj[P_SITE_POS1] = getSiteData($$aln[S_RNAME], $proj[P_SITE_INDEX1]);
    @proj;
}
sub getSiteData {
    my ($chrom, $siteIndex1) = @_;
    my $siteDataOffset = ${$chromData{$chrom}}[siteDataOffset_] + ($siteIndex1 - 1) * BYTES_PER_SITE_DATA;
    seek($siteDataH, $siteDataOffset, SEEK_SET) or die "could not seek in siteData lookup: $!\n";
    if(!read($siteDataH, $siteDataRaw, BYTES_PER_SITE_DATA)){
        # print STDERR "$error: could not read from siteIndex lookup: $chrom, $adjPos1: $!\n";
        return 0;
    }
    my ($sitePos1) = unpack(SITE_DATA_PACKING, $siteDataRaw);
    $sitePos1; # discard fragment sizes here
}

# clean up
sub finishMatchSites {
    close $siteIndexH;
    close $siteDataH;
}

1;

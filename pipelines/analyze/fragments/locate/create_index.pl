use strict;
use warnings;

# pack index files for fast retrieval of the RE sites closest to alignment endpoints
# three output files are created:
#   SITE_CHROM_DATA_FILE      = table of chrom offsets for the next two lookup files
#       one row per chrom with columns: chrom,chromIndex,nSites,chromSize,closestSiteOffset,siteDataOffset
#   CLOSEST_SITE_BINARY_FILE  = one value per refPos1 = signed closest siteIndex1 (+ means refPos1 >= closestSitePos1)
#   SITE_DATA_BINARY_FILE     = one value set per RE filtering site = sitePos1 on chrom, two flanking fragment sizes
# the use of site indices along a chrom as an intervening step to sitePos1 allows:
#   assignment to numbered RE fragments 
#   calculation of distances to closest sites to enforce site tolerance
#   quickly jumping from one site to the next by index, not position
#   use of smaller values to store site information per read
# finally, the filtering sites stream is also written to STDOUT for tabix indexing for app

# initialize reporting
our $action = "create_index";
my ($nInputFiles) = (0) x 10;

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);

# environment variables
fillEnvVar(\our $GENOME_FASTA,               'GENOME_FASTA');
fillEnvVar(\our $EXPECTING_ENDPOINT_RE_SITES,'EXPECTING_ENDPOINT_RE_SITES');
# sample-level input and output files
fillEnvVar(\our $FILTERING_SITES_FILE,      'FILTERING_SITES_FILE');
fillEnvVar(\our $SITE_CHROM_DATA_FILE,      'SITE_CHROM_DATA_FILE');
fillEnvVar(\our $CLOSEST_SITE_BINARY_FILE,  'CLOSEST_SITE_BINARY_FILE');
fillEnvVar(\our $SITE_DATA_BINARY_FILE,     'SITE_DATA_BINARY_FILE');
# genome-level input and output files
fillEnvVar(\our $GENOME_FILTERING_SITES_FILE,      'GENOME_FILTERING_SITES_FILE');
fillEnvVar(\our $GENOME_SITE_CHROM_DATA_FILE,      'GENOME_SITE_CHROM_DATA_FILE');
fillEnvVar(\our $GENOME_CLOSEST_SITE_BINARY_FILE,  'GENOME_CLOSEST_SITE_BINARY_FILE');
fillEnvVar(\our $GENOME_SITE_DATA_BINARY_FILE,     'GENOME_SITE_DATA_BINARY_FILE');
my $expectingEndpointReSites = $EXPECTING_ENDPOINT_RE_SITES eq "TRUE" ? 1 : 0;
unless($expectingEndpointReSites){ # $REJECTING_JUNCTION_RE_SITES is always true if this script was called
    $FILTERING_SITES_FILE     = $GENOME_FILTERING_SITES_FILE;
    $SITE_CHROM_DATA_FILE     = $GENOME_SITE_CHROM_DATA_FILE;
    $CLOSEST_SITE_BINARY_FILE = $GENOME_CLOSEST_SITE_BINARY_FILE;
    $SITE_DATA_BINARY_FILE    = $GENOME_SITE_DATA_BINARY_FILE;
}

# initialize the genome
use vars qw(%chromIndex %revChromIndex @canonicalChroms);
setCanonicalChroms();

# constants
use constant {
    NAME        => 0, # FAI fields
    LENGTH      => 1,
    OFFSET      => 2,
    LINEBASES   => 3, # 1-based
    LINEWIDTH   => 4,
    #-----------------
    CHROM       => 0, # filtering sites table columns
    SITE_POS1   => 1,
    # IN_SILICO   => 2,
    # N_OBSERVED  => 3,
    #-----------------
    CLOSEST_SITE_PACKING   => "l", # signed integer (siteIndex1)
    BYTES_PER_CLOSEST_SITE => 4,
    SITE_DATA_PACKING      => "LSS", # unsigned integer (sitePos1), unsigned shorts (flanking frag sizes)
    BYTES_PER_SITE_DATA    => 8,
    ITEMS_PER_PACK_CALL    => 10000,
    MAX_SHORT_VALUE        => 65535,
    #-------------
    LEFTWARD  => -1,
    RIGHTWARD =>  1,
};

# load chromosome sizes (all positions will be filled in CLOSEST_SITE_BINARY_FILE)
print STDERR "  loading chrom sizes\n";
my (%chromSizes);
open my $faiH, "<", "$GENOME_FASTA.fai" or die "file not found: $GENOME_FASTA.fai\n";
while(my $line = <$faiH>){
    chomp $line;
    my @chrom = split("\t", $line);
    $chromSizes{$chrom[NAME]} = $chrom[LENGTH];
}
close $faiH; 

# load the RE filtering sites
print STDERR "  loading filtering sites\n";
my (%sitePos1, %nSites, $chr, $chrIdx1);
open my $sitesH, "-|", "zcat $FILTERING_SITES_FILE" or die "could not open: $!\n";
if($expectingEndpointReSites){
    my $header = <$sitesH>; # FILTERING_SITES_FILE has a header line, GENOME_FILTERING_SITES_FILE==GENOME_SITES_GZ does not
}
while (my $line = <$sitesH>){ # chrom,sitePos1[,inSilico,nObserved]
    chomp $line;
    my @site = split("\t", $line);
    ($chr, $chrIdx1) = ($site[CHROM], $chromIndex{$site[CHROM]});
    push @{$sitePos1{$chr || "NA"}}, $site[SITE_POS1];
    $chrIdx1 or next;
    $site[CHROM] = $chrIdx1;
    print join("\t", @site), "\n";
}
close $sitesH;
foreach my $chrom(keys %sitePos1){
    $nSites{$chrom} = scalar(@{$sitePos1{$chrom}});
}

# write the tabular index file
print STDERR "  writing chrom data with offsets\n";
open my $indexH, ">", $SITE_CHROM_DATA_FILE or die "could not open: $!";
my ($closestSiteOffset, $siteDataOffset) = (0, 0);
print $indexH join("\t", qw(chrom chromIndex nSites chromSize closestSiteOffset siteDataOffset)), "\n";
foreach my $chrom(@canonicalChroms){
    $nSites{$chrom} or $nSites{$chrom} = 0;
    print $indexH join("\t",
        $chrom, 
        $chromIndex{$chrom},
        $nSites{$chrom},
        $chromSizes{$chrom},
        $closestSiteOffset, # this genome-scaled value requires 64-bit integers to load (all others are 32-bit)
        $siteDataOffset
    ), "\n";
    $closestSiteOffset +=  $chromSizes{$chrom} * BYTES_PER_CLOSEST_SITE;
    $siteDataOffset    +=  $nSites{$chrom}     * BYTES_PER_SITE_DATA;
}
close $indexH;

# write the binary lookup files
open my $siteIndexH, ">:raw", $CLOSEST_SITE_BINARY_FILE or die "could not open file: $!";
open my $siteDataH,  ">:raw", $SITE_DATA_BINARY_FILE    or die "could not open file: $!";

# commit one chromosome at a time
print STDERR "  writing binary lookup files\n";
my ($sitePos1, $maxSiteI);
foreach my $chrom(@canonicalChroms){
    print STDERR join("    ", "", $chrom, commify($nSites{$chrom})." sites", commify($chromSizes{$chrom})." bases"), "\n";   
    $sitePos1 = $sitePos1{$chrom};
    ($sitePos1 and $nSites{$chrom}) or next;
    $maxSiteI = $nSites{$chrom} - 1;

    # commit all chrom refPos1 before the first sitePos1 as negative index
    # site indices are 1-referenced to always be able to take a sign
    my $i0 = 0;
    my $i1 = $i0 + 1;
    my $n = $$sitePos1[$i0] - 1;
    writeClosestSites(-$i1, $n);

    # process all sites on the chrom
    my $lastCommitedPos1;
    foreach my $i0(0..$#{$sitePos1}){
        my $i1 = $i0 + 1;

        # commit half the distance from the previous site as negative index
        if($i0 > 0){
            my $n = $$sitePos1[$i0] - $lastCommitedPos1 - 1;
            writeClosestSites(-$i1, $n);
        }

        # commit sitePos1 itself 
        writeClosestSites($i1, 1); # the lookup from sitePos1 to siteIndex1 of the closest site
        print $siteDataH pack(
            SITE_DATA_PACKING,               # the reverse lookup from siteIndex1 to ...
            $$sitePos1[$i0],                 # ... sitePos1
            getFragmentSize($i0,  LEFTWARD), #  leftward fragment size
            getFragmentSize($i0, RIGHTWARD)  # rightward fragment size
        );

        # commit half the distance to the next site as positive index
        if($i0 < $#{$sitePos1}){
            my $n = int(($$sitePos1[$i0 + 1] - $$sitePos1[$i0] - 1) / 2);
            writeClosestSites($i1, $n);
            $lastCommitedPos1 = $$sitePos1[$i0] + $n;
        }
    }

    # commit all chrom refPos1 after the last sitePos1 as positive index
    $i0 = $#{$sitePos1};
    $i1 = $i0 + 1;
    $n = $chromSizes{$chrom} - $$sitePos1[$i0];
    writeClosestSites($i1, $n)
}
print STDERR "  finishing up\n";

# clean up
close $siteIndexH;
close $siteDataH;

# iterate pack just in case there are very large gaps
sub writeClosestSites {
    my ($value, $n) = @_;
    while($n >= ITEMS_PER_PACK_CALL){
        print $siteIndexH pack(CLOSEST_SITE_PACKING x ITEMS_PER_PACK_CALL, ($value) x ITEMS_PER_PACK_CALL);
        $n -= ITEMS_PER_PACK_CALL;
    }
    $n or return;
    print $siteIndexH pack(CLOSEST_SITE_PACKING x $n, ($value) x $n);
}

# calculate the distance from an index sitePos1 to the next sitePos1
# work in either the leftward or rightward direction from the index sitePos1 as requested
sub getFragmentSize {
    my ($i0, $direction) = @_;
    $direction ==  LEFTWARD and $i0 == 0         and return 0; # telomeric fragments are not used
    $direction == RIGHTWARD and $i0 == $maxSiteI and return 0;
    my $j0 = $i0 + $direction;
    min(
        MAX_SHORT_VALUE, # thus, support RE fragments up to 65K on long-read platforms
        abs(
            $$sitePos1[$i0] - 
            $$sitePos1[$j0]
        )
    );
}

use strict;
use warnings;

# action:
#   collect all reference genome RE sites for selected blunt restriction enzymes
#   handle degenerate sites using regular expressions in the RE sites table
#   handle non-palindromic sites by paired site entries in the RE sites table
# input: 
#   GENOME_FASTA on STDIN
# output:
#   one compressed RE site position file per RE

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/genome/$_.pl" } qw(chroms);

# environment variables
fillEnvVar(\our $GENOME_REMAPS_DIR, 'GENOME_REMAPS_DIR');
fillEnvVar(\our $GENOME,            'GENOME');
fillEnvVar(\our $BLUNT_RE_TABLE,    'BLUNT_RE_TABLE');
fillEnvVar(\our $MIN_PRIORITY_LEVEL,'MIN_PRIORITY_LEVEL');

# initialize the genome
use vars qw(%chromIndex);
setCanonicalChroms();

# initialize the REs
# enzyme  strand  cut_site regex   offset  CpG_priority
# EcoRV   0       GATATC   GATATC  3       4  
# thus:
#            *
#  EcoRV GAT^ATC
#        CTA^TAG
#          *
print STDERR "using restriction enyzmes:\n";
my (@reSites, @reKeys, @reOffsets, @reRegex);
open my $inH, "<", $BLUNT_RE_TABLE or die "could not open: $!\n";
my $header = <$inH>; # enzyme,strand,cut_site,regex,offset,CpG_priority,...
while (my $line = <$inH>){
    my ($enzyme, $strand, $cut_site, $regex, $offset, $priority) = split(",", $line);
    $enzyme or next;
    $priority < $MIN_PRIORITY_LEVEL and next;
    $enzyme =~ s/\s//g;
    $cut_site =~ s/\s//g;
    $cut_site = uc($cut_site);
    push @reSites, $cut_site;
    push @reKeys, "$enzyme\_$cut_site"; # non-palindromic cut sites generate two files, one per strand
    push @reOffsets, $offset;
    push @reRegex, qr/$regex/;
    print STDERR "  $enzyme\n";
}
close $inH;
my @reIs = 0..$#reSites;

# initialize the process
my ($prevLen, $prevCumLength, $chrom, $prevLine) = (0, 0);
my ($genomeSize, %siteCounts);
my @fileHs = map {
    my $file = "$GENOME_REMAPS_DIR/$GENOME.digest.$_.txt.gz";
    open my $fileH, "|-", "gzip -c > $file" or die "could not open stream: $file\n";
    $fileH;
} @reKeys;

# handle site finding
sub findSites {
    my ($line) = @_;
    my $query = uc("$prevLine$line");
    foreach my $reI(@reIs){
        my $regex = $reRegex[$reI];
        while($query =~ m/$regex/g){
            my $queryPos = $-[0];
            $queryPos < $prevLen or last;

            # the 1-indexed first base AFTER the cleaved phosphodiester bond on the top strand of the genome
            my $sitePos1 = $prevCumLength + $queryPos + $reOffsets[$reI] + 1;
            my $fH = $fileHs[$reI];
            print $fH join("\t", $chrom, $sitePos1), "\n";
            $siteCounts{$reKeys[$reI]}++;
        }
    }
}

# process the genome fasta file one line at a time
print STDERR "performing in silico digestions:\n";  
while (my $line = <STDIN>){
    chomp $line;
    if($line =~ m/^>(\S+)/){
        $chrom and $prevLine and findSites("");
        $chrom = $1;
        $chromIndex{$chrom} or $chrom = ""; # don't include non-canonical chroms
        $prevLine = "";
        $prevLen = 0;
        $prevCumLength = 0;
        $chrom and print STDERR "  $chrom\n";
    } elsif($chrom) {
        $prevLine and findSites($line);
        $prevCumLength += $prevLen;
        $prevLine = $line;
        $prevLen = length($prevLine);
        $genomeSize += $prevLen;
    };
}
$chrom and $prevLine and findSites("");
foreach my $fileH(@fileHs){
    close $fileH;
}

# print the final tallies
print STDERR "total site counts in ".commify($genomeSize)." genome bp\n";
foreach my $reKey(@reKeys){
    print STDERR join("\t", "", $reKey, commify($siteCounts{$reKey} || 0)), "\n";
}

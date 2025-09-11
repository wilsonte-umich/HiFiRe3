use strict;
use warnings;

#----------------------------------------------------------
# genome manipulations
#----------------------------------------------------------

# working variables
our (@canonicalChroms, %chromIndex, %revChromIndex);

# parsers for restricting work to properly ordered canonical chromosomes
sub setCanonicalChroms { 
    my @chroms = split(/\s+/, $ENV{GENOME_CHROMS}); # all placed chromosome sequences including chrM and chrEBV if present (but not chrXX_XX)
    if($ENV{USE_ALL_CHROMS} or $ENV{IS_COMPOSITE_GENOME}){
        @canonicalChroms = @chroms; # chroms used in fai indexed order
    } else {
        my %canonicalChroms = map { $_ => 1 } @chroms;
        sub getPushValue{
            my ($chr, $canonicalChroms) = @_;
            my $chrom = "chr$chr";
            $$canonicalChroms{$chrom} ? $chrom : ();
        } 
        my $isRoman = (uc($ENV{GENOME_CHROMS}) =~ m/CHRI/);
        @canonicalChroms = $isRoman ? 
            map { getPushValue($_, \%canonicalChroms) } qw(I II III IV V VI VII VIII IX X 
                                                            XI XII XIII XIV XV XVI XVII XVIII XIX XX 
                                                            XXI XXII XXIII XXIV XXV XXVI XXVII XXVIII XXIX XXX) :
            map { getPushValue($_, \%canonicalChroms) } (1..90,'X','Y','M','EBV');
    }
    %chromIndex    = map { $canonicalChroms[$_] => $_ + 1 } 0..$#canonicalChroms; # 1-referenced chrom indices, i.e., chr3 => 3
    %revChromIndex = map { $_ + 1 => $canonicalChroms[$_] } 0..$#canonicalChroms;    
    $chromIndex{'*'} = 99; # special handling of unmapped reads
    $revChromIndex{99} = '*'; 
}
sub getNuclearChroms {
    map {
        (uc($_) =~ m/CHRM/ or uc($_) =~ m/CHREBV/)? () : $_;
    } @canonicalChroms;
}
sub getChromSizes {
    my ($faiFile) = @_;
    my %chromSizes;
    open my $faiH, "<", $faiFile or die "could not open: $faiFile: $!\n";
    while(my $line = <$faiH>){
        chomp $line;
        my ($chrom, $size) = split("\t", $line);
        $chromSizes{$chrom} = $size;
    }
    close $faiH;
    return %chromSizes;
}
sub getChromIndexSizes {
    my ($faiFile) = @_;
    my @chromSizes;
    open my $faiH, "<", $faiFile or die "could not open: $faiFile: $!\n";
    while(my $line = <$faiH>){
        chomp $line;
        my ($chrom, $size) = split("\t", $line);
        $chromIndex{$chrom} and $chromSizes[$chromIndex{$chrom}] = $size;
    }
    close $faiH;
    return @chromSizes;
}
sub writeChromsFile {
    my ($chromsFile, $genomeFasta) = @_;
    @canonicalChroms or setCanonicalChroms();
    my @chromSizes = getChromIndexSizes("$genomeFasta.fai");
    open my $chrH, ">", $chromsFile or die "could not open: $chromsFile: $!\n";
    foreach my $chrom(keys %chromIndex){
        my $chromIndex1 = $chromIndex{$chrom};
        print $chrH join("\t", $chrom, $chromIndex1, $chromSizes[$chromIndex1] || "NA"), "\n";
    }
    close $chrH;
}

1;

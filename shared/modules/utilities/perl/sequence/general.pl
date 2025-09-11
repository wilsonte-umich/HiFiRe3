use strict;
use warnings;

#----------------------------------------------------------
# general sequence manipulations
#----------------------------------------------------------

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
    READ_SV_PENALTY_FRACTION => 0.2, # SV penalty as a function of read length
};

# regular expressions to process read clips
our $leftClip_  = qr/^(\d+)S/;
our $rightClip_ = qr/(\d+)S$/;

# conversion from SAM encoding to PAF encoding
sub getQueryStart0 {
    my ($flag, $cigar) = @_;
    ($flag & _UNMAPPED) and return 0;
    if($flag & _REVERSE){
        $cigar =~ m/$rightClip_/ ? $1 : 0;
    } else {
        $cigar =~ m/$leftClip_/  ? $1 : 0;
    }
}
sub getQueryEnd1 {
    my ($flag, $cigar, $qLen) = @_;
    ($flag & _UNMAPPED) and return $qLen;
    if($flag & _REVERSE){
        $qLen - ($cigar =~ m/$leftClip_/  ? $1 : 0);
    } else {
        $qLen - ($cigar =~ m/$rightClip_/ ? $1 : 0);
    }
}

# get SAM clips lengths from input CIGAR string
sub getLeftClip{
    $_[0] =~ m/$leftClip_/  ? $1 : 0;
}
sub getRightClip{
    $_[0] =~ m/$rightClip_/ ? $1 : 0;
}

# get rightmost mapped read pos1 in reference genome from POS and CIGAR
sub getEnd { 
    my ($pos1, $cigar) = @_;
    $cigar =~ s/\d+S//g;
    my $end = $pos1 - 1;
    while ($cigar =~ (m/(\d+)(\w)/g)) {
        $2 eq "I" or $end += $1;
    }
    return $end;
}
sub getAlignedSize {
    my ($cigar) = @_;
    $cigar =~ s/\d+[HSD]//g;
    my $size = 0;
    while ($cigar =~ (m/(\d+)(\w)/g)) {
        $size += $1;
    }
    return $size;
}
sub getRefSpan {
    my ($cigar) = @_;
    $cigar =~ s/\d+[HSI]//g;
    my $size = 0;
    while ($cigar =~ (m/(\d+)(\w)/g)) {
        $size += $1;
    }
    return $size;
}

# reverse complement reads, by reference ...
sub rc { 
    my ($seq) = @_;
    $$seq = reverse $$seq;
    $$seq =~ tr/ATGC/TACG/;
}
#... and by value
sub getRc {
    my ($seq) = @_;
    rc(\$seq);
    $seq;
}

# get the average per-base Phred QUAL for a single read or segment of a read
sub getAvgQual {
    ($_[0] and length($_[0])) or return 0;
    my $sum = 0;
    map{ $sum += ord($_) } split("", $_[0]);
    $sum / length($_[0]) - 33;
}
sub getAvgQual_round5 {
    roundCount(getAvgQual($_[0]), 0.2);
}

# parse a CIGAR string to match a query SEQ to its reference
sub getQryOnRef {
    my ($qry, $cigar, $noClip) = @_;
    my @qry = split("", $qry);
    my @qryOnRef;
    my $qryI = 0;
    my $nDeleted = 0;
    my $nInserted = 0;
    my @insI;
    if($noClip){
        $cigar =~ s/^(\d+)[S|H]//g and @qry = @qry[$1..$#qry];
        $cigar =~ s/(\d+)[S|H]$//g and @qry = @qry[0..($#qry - $1)];
    }
    while ($cigar =~ (m/(\d+)(\w)/g)) { 
        my ($size, $operation) = ($1, $2);
        if($operation eq 'D'){
            push @qryOnRef, (("-") x $size);            
            $nDeleted += $size;
        } elsif($operation eq 'I'){
            push @insI, $qryI - 1 + $nDeleted - $nInserted;
            $nInserted += $size;
            $qryI += $size; 
        } else {
            push @qryOnRef, @qry[$qryI..($qryI + $size - 1)];
            $qryI += $size;
        } 
    }
    foreach my $i(@insI){ $qryOnRef[$i] = "+" } # mark the position to the left of each novel insertion
    return \@qryOnRef;
}

# analyze a length distribution
sub getN50 {
    my (@sizes) = @_;
    my $sumSizes = 0;
    map { $sumSizes += $_ } @sizes;
    my $halfSumSizes = $sumSizes / 2;
    my $runningSum = 0;
    foreach my $size(sort {$b <=> $a} @sizes){
        $runningSum += $size;
        $runningSum >= $halfSumSizes and return $size;
    }
}

# aggregate alignment scores over all segments of a read (e.g., one half of a paired read sequence)
# invoke an SV junction/split-read penalty to help suppress reference bias that leads to SV artifacts
#----------------------------------------------------------------------------------
# minimap2 in short-read mode has a match score of 2 and a mismatch penalty of 8, for a net penalty of 10
# so a 150 bp read with zero and one mismatches have scores of 300 and 290, respectively
# a READ_SV_PENALTY_FRACTION of 0.2 means a 150 bp read SV would offset 3 mismatches (150 * 0.2 = 30 = 3 * 10)
#----------------------------------------------------------------------------------
# minimap2 in standard long-read modes has a match score of 2 and a mismatch penalty of 4, for a net penalty of 6
sub getCompositeReadAlnScore {
    my ($alns) = @_;
    my $scoreSum = 0;
    foreach my $aln(@$alns){
        $scoreSum += (($$aln[TAGS] and $$aln[TAGS] =~ m/AS:i:(\S+)/) ? $1 : 0);
    }
    $scoreSum - (@$alns - 1) * READ_SV_PENALTY_FRACTION * length($$alns[0][SEQ]);
}

# get the minimum max-insert-size from config and options
sub getMaxInsertSize {
    $ENV{PLATFORM_MAX_INSERT_SIZE} or die "missing or zero value for PLATFORM_MAX_INSERT_SIZE\n";
    $ENV{MAX_INSERT_SIZE}          or die "missing or zero value for MAX_INSERT_SIZE\n";
    min($ENV{PLATFORM_MAX_INSERT_SIZE}, $ENV{MAX_INSERT_SIZE});
}

1;

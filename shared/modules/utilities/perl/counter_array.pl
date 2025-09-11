use strict;
use warnings;

# memory efficient counter array for a large number of positions, e.g., bases on a chromosome
# increments and returns counts for each array index, e.g., for each base position on a chromosome
# each index is a 16-bit integer that uses 2 bytes of memory and can store counts up to 65,535

use constant {
    BYTES_PER_INTEGER => 2,   # for 16-bit unsigned int
    N_POS_IN_CHUNK    => 1e6, # 1M positions at a time
    N_BYTES_IN_CHUNK  => 2e6,
};

sub counter_array_init {
    my ($arrayLength) = @_;
    my $data = '';
    open my $fh, '>', \$data or die "counter_array_init failed: $!";
    my $zeros = "\0" x (N_POS_IN_CHUNK * BYTES_PER_INTEGER);
    my $nPosRemaining = $arrayLength;
    while($nPosRemaining > 0) {
        print $fh $zeros; # will overrun the required buffer by up to 1M positions, 2MB memory
        $nPosRemaining -= N_POS_IN_CHUNK;
    }
    close $fh;
    open $fh, '+<', \$data or die "counter_array_init failed: $!";
    {
        fh => $fh,
        arrayLength => $arrayLength,
        dataRef => \$data
    }
}
sub counter_array_init_from_file {
    my ($filename) = @_;
    my $data = '';
    open my $fh, '>', \$data or die "counter_array_init_from_file open 1 failed: $!";
    open my $inH, '<', $filename or die "counter_array_init_from_file open 2 failed: $!";
    my ($buffer, $arrayLength);
    while (sysread($inH, $buffer, N_BYTES_IN_CHUNK)) {
        print $fh $buffer;
        $arrayLength += length($buffer) / BYTES_PER_INTEGER;
    }
    close $fh;
    close $inH;
    open $fh, '+<', \$data or die "counter_array_init failed: $!";
    {
        fh => $fh,
        arrayLength => $arrayLength,
        dataRef => \$data
    }
}

sub counter_array_get {
    my ($counter, $index0) = @_;
    $index0 < $$counter{arrayLength} or return 0;
    my $fh = $$counter{fh};
    seek($fh, $index0 * BYTES_PER_INTEGER, 0) or die "counter_array_get seek failed: $!";
    my $value;
    read($fh, $value, BYTES_PER_INTEGER) or die "counter_array_get read failed: index $index0: $!";
    unpack('S', $value);
}
sub counter_array_set {
    my ($counter, $index0, $value) = @_;
    $index0 < $$counter{arrayLength} or return;
    my $fh = $$counter{fh};
    seek($fh, $index0 * BYTES_PER_INTEGER, 0) or die "counter_array_set seek failed: $!";
    print $fh pack('S', $value);
}
sub counter_array_set_range {
    my ($counter, $startIndex0, $length, $value) = @_;
    $startIndex0 + $length <= $$counter{arrayLength} or return;
    my $fh = $$counter{fh};
    seek($fh, $startIndex0 * BYTES_PER_INTEGER, 0) or die "counter_array_set seek failed: $!";
    print $fh pack('S', $value) x $length; # sets the same value into a range of positions
}
sub counter_array_increment {
    my ($counter, $index0, $count) = @_;
    $index0 < $$counter{arrayLength} or return;
    my $value = counter_array_get($counter, $index0);
    counter_array_set($counter, $index0, $value + $count);
}

sub counter_array_close {
    my ($counter) = @_;
    my $fh = $$counter{fh};
    close $fh;
}
sub counter_array_write_to_file {
    my ($counter, $filename) = @_;
    open my $outH, '>:raw', $filename or die "counter_array_write_to_file failed: $!";
    my $dataRef = $$counter{dataRef};
    print $outH $$dataRef; # prints packed data in binary format for long-term storage
    close $outH;
}
sub counter_array_release {
    my ($counter, $close) = @_;
    $close and counter_array_close($counter);
    my $dataRef = $$counter{dataRef};
    $$dataRef = ''; 
}

# debugging support
sub counter_array_size {
    my ($counter, $debug) = @_;
    my $dataRef = $$counter{dataRef};
    my $dataRefLen = length($$dataRef) / BYTES_PER_INTEGER;
    print STDERR "counter_array_size: $debug: $$counter{arrayLength}/$dataRefLen\n";
    foreach my $offsetMultiplier(0..5){
        my $sample = substr($$dataRef, 10e6 * $offsetMultiplier, 50 * BYTES_PER_INTEGER);
        $sample = join(",", unpack('S*', $sample));
        print STDERR "counter_array_size: $debug: $sample\n";
    }
    my @x;
    for my $i(1..50){
        push @x, counter_array_get($counter, int(rand($$counter{arrayLength} - 1)));
    }
    print STDERR join(",", @x), "\n";
}

1;

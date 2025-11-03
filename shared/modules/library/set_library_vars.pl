use strict;
use warnings;

# action (in order):
#   read the default variables from the config file
#   override defaults with non-null user options
#   cascade to set derived variables
#   print for sourcing by the calling bash script

# set working variables
my %envVars;

# set the working platform
my @platforms = qw(ONT PacBio);
my $platformI;
foreach my $i(0..$#platforms){
    if($platforms[$i] eq $ENV{SEQUENCING_PLATFORM}){
        $platformI = $i;
        last;
    }
}
unless(defined $platformI){
    print "export GET_CONFIGURATION_ERROR='unrecognized --sequencing-platform: $ENV{SEQUENCING_PLATFORM}'\n";
    exit 1;
}

# read defaults
# override defaults with user job options
my $inH;
unless(open $inH, "<", "$ENV{MODULES_DIR}/library/libraries.csv"){
    print "export GET_CONFIGURATION_ERROR='could not open libraries.csv: $!'\n";
    exit 1;
}
my $line = <$inH>; # discard header
while (my $line = <$inH>){
    my ($category, $envVar, @values) = split(",", $line); # 0 or NA values are ignored and not used as overrides
    $values[$platformI] =~ s/;/,/g;
    $envVars{$envVar} = 
        ($ENV{$envVar} and $ENV{$envVar} ne "NA" and $ENV{$envVar} ne "null") ? 
        $ENV{$envVar} :
        $values[$platformI];
}
close $inH;

# set derived environment variables
# nothing...

# print for sourcing by calling script
foreach my $envVar(keys %envVars){
    $envVar or next;
    print "export $envVar='$envVars{$envVar}'\n";
}

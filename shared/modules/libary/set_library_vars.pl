use strict;
use warnings;

# action (in order):
#   read the platform default variables from the default config file
#   override defatuls with non-null user options
#   cascade to set derived variables
#   print for sourcing by the calling bash script

# set working variables
my %envVars;

# set the working platform
my @platforms = qw(IlluminaPE AvitiPE AvitiSE Ultima ONT PacBio);
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
unless(open $inH, "<", "$ENV{MODULES_DIR}/agFree/platforms.csv"){
    print "export GET_CONFIGURATION_ERROR='could not open platforms.csv: $!'\n";
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
$envVars{FASTP_INTERLEAVED_OPTIONS} = 
    $envVars{READ_PAIR_TYPE} eq "paired" ? 
    "--interleaved_in":
    "";
$envVars{FASTP_MERGE_OPTIONS} = 
    (
        $envVars{READ_PAIR_TYPE} eq "paired" and 
        $envVars{MERGE_READS}
    ) ? 
    "--merge --include_unmerged":
    "";
$envVars{FASTP_MERGE_NOTICE} = $envVars{FASTP_MERGE_OPTIONS} ? "with pre-alignment merging" : "";
$envVars{FASTP_POLY_X_OPTIONS} = 
    $envVars{TRIM_POLY_X} ?
        "--trim_poly_x --poly_x_min_len $envVars{POLY_X_MIN_LEN}" : 
        "";

# print for sourcing by calling script
foreach my $envVar(keys %envVars){
    $envVar or next;
    print "export $envVar='$envVars{$envVar}'\n";
}

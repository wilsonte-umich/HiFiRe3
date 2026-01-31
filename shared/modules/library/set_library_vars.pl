use strict;
use warnings;

# action (in order):
#   validate the values for options --sequencing-platform and --library-type
#   read the platform and library type default variables from the config files
#   override defaults with non-null user options
#   check that the library type is supported for the selected platform
#   cascade to set derived variables
#   print for sourcing by the calling bash script

# set working variables
my %envVars;

# set and check the requested platform and library type
sub getOptionI {
    my ($option, $envKey, @values) = @_;
    foreach my $i(0..$#values){
        if($values[$i] eq $ENV{$envKey}){
            return $i;
        }
    }
    print "export GET_CONFIGURATION_ERROR='unrecognized $option: $ENV{$envKey}'\n";
    exit 1;
}
my $platformI = getOptionI(
    "--sequencing-platform",
    "SEQUENCING_PLATFORM", 
    qw(Illumina_2x150 Aviti_2x150 Aviti_1x300 Ultima ONT PacBio)
);
my $libTypeI = getOptionI(
    "--library-type",
    "LIBRARY_TYPE", 
    qw(Nextera TruSeq Elevate Ultima Ligation Rapid HiFi)
);

# read default option values for the selected platform and library type
# override defaults with user job options
sub applyDefaults {
    my ($fileName, $index) = @_;
    my $inH;
    unless(open $inH, "<", "$ENV{MODULES_DIR}/library/$fileName.csv"){
        print "export GET_CONFIGURATION_ERROR='could not open $fileName.csv: $!'\n";
        exit 1;
    }
    my $line = <$inH>; # discard header
    while (my $line = <$inH>){
        my ($category, $envVar, @values) = split(",", $line); # 0 or NA values are ignored and not used as overrides
        $values[$index] =~ s/;/,/g;
        $envVars{$envVar} = 
            ($ENV{$envVar} and $ENV{$envVar} ne "NA" and $ENV{$envVar} ne "null") ? 
            $ENV{$envVar} :
            $values[$index];
    }
    close $inH;
}
applyDefaults("platforms", $platformI);
applyDefaults("libraries", $libTypeI);

# check that the library type is compatible with the platform
unless($envVars{SUPPORTED_PLATFORMS} =~ /\|$ENV{SEQUENCING_PLATFORM}\|/){
    print "export GET_CONFIGURATION_ERROR='library type $ENV{LIBRARY_TYPE} is not supported for platform $ENV{SEQUENCING_PLATFORM}'\n";
    exit 1;
}

# set derived environment variables
$envVars{RUN_PREALIGNMENT_FASTP} = 
    ($envVars{READ_PAIR_TYPE} eq "paired" or $envVars{TRIM_POLY_X} or $envVars{TRIM_ADAPTERS}) ? 
    "TRUE" : 
    "";
$envVars{FASTP_PAIRED_END_OPTIONS} = 
    $envVars{READ_PAIR_TYPE} eq "paired" ? 
    "--interleaved_in --merge --include_unmerged":
    "";
$envVars{FASTP_POLY_X_OPTIONS} = 
    $envVars{TRIM_POLY_X} ?
    "--trim_poly_x --poly_x_min_len $envVars{POLY_X_MIN_LEN}" : 
    "";
$envVars{FASTP_ADAPTER_OPTIONS} = 
    $envVars{TRIM_ADAPTERS} ? join(" ", 
        (
            $envVars{ADAPTER_SEQUENCE} ? 
            "--adapter_sequence $envVars{ADAPTER_SEQUENCE}" : 
            ""
        ),
        (
            $envVars{ADAPTER_SEQUENCE_R2} ? 
            "--adapter_sequence_r2 $envVars{ADAPTER_SEQUENCE_R2}" : 
            ""
        )
    ) : "--disable_adapter_trimming";
$envVars{DEDUPLICATE_READS} = (
    $envVars{DEDUPLICATE_PLATFORM} or 
    $envVars{DEDUPLICATE_LIBRARY}
) ? "TRUE" : "";

# print for sourcing by calling script
foreach my $envVar(keys %envVars){
    $envVar or next;
    print "export $envVar='$envVars{$envVar}'\n";
}

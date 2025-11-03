use strict;
use warnings;

# actions:
#   merge event metadata from junction_sources into junctions
# input:
#   $SV_UNIQUE_JUNCTIONS_FILE
#   $SV_JUNCTION_SOURCES_FILE
# output:
#   unique junctions with merged source metadata on STDOUT

# initialize reporting
our $script = "merge_events";
our $error  = "$script error";

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);

# environment variables
fillEnvVar(\our $SV_UNIQUE_JUNCTIONS_FILE, 'SV_UNIQUE_JUNCTIONS_FILE');
fillEnvVar(\our $SV_JUNCTION_SOURCES_FILE, 'SV_JUNCTION_SOURCES_FILE');

# constants
use constant {
    OJXN_CHROM_INDEX1_1     => 0, # leader columns shared and sorted the same in SV_UNIQUE_JUNCTIONS_FILE and SV_JUNCTION_SOURCES_FILE
    OJXN_REF_POS1_1         => 1,
    OJXN_STRAND_INDEX0_1    => 2,
    OJXN_CHROM_INDEX1_2     => 3,
    OJXN_REF_POS1_2         => 4,
    OJXN_STRAND_INDEX0_2    => 5,
    OJXN_JXN_TYPE           => 6,
    #-------------
    UJXN_N_OBSERVED         => 7, # continuation columns in SV_UNIQUE_JUNCTIONS_FILE
    UJXN_ALN_OFFSET         => 8,
    UJXN_JXN_BASES          => 9,
    UJXN_PATHS              => 10,
    UJXN_N_PATH_JUNCTIONS   => 11,
    UJXN_JSRC_MAPQ          => 12, # added by this script, a single junction_sources value for the selected, most frequent junction properties
    UJXN_JSRC_DE_TAG        => 13,
    UJXN_JSRC_SITE_DIST     => 14,
    UJXN_JSRC_ALN_FAILURE_FLAG => 15,
    UJXN_JSRC_JXN_FAILURE_FLAG => 16,
    UJXN_JSRC_QNAMES           => 17, # added by this script, concatenated lists across all junction_sources
    UJXN_JSRC_INSERT_SIZES     => 18,
    UJXN_JSRC_IS_ALLOWED_SIZES => 19,
    UJXN_JSRC_MIN_STEM_LENGTHS => 20,
    UJXN_JSRC_PASSED_STEMS     => 21,
    UJXN_JSRC_SEQS            => 22,
    UJXN_JSRC_QUALS           => 23,
    UJXN_JSRC_CIGARS          => 24,
    UJXN_JSRC_ORIENTATIONS    => 25,
    UJXN_JSRC_CHANNELS        => 26, # concatenated lists that are dropped later after used for duplicate purging
    UJXN_JSRC_OUTER_ENDPOINTS => 27,
    #-------------
    JSRC_ALN_OFFSET         => 7, # continuation columns in SV_JUNCTION_SOURCES_FILE
    JSRC_JXN_BASES          => 8,
    JSRC_ORIENTATIONS       => 9,
    JSRC_QNAMES             => 10,
    JSRC_CHANNELS           => 11,
    JSRC_INSERT_SIZES       => 12,
    JSRC_IS_ALLOWED_SIZES   => 13,
    JSRC_OUTER_ENDPOINTS    => 14,
    JSRC_MIN_STEM_LENGTHS   => 15,
    JSRC_PASSED_STEMS       => 16,
    JSRC_MAPQ               => 17,
    JSRC_DE_TAG             => 18,
    JSRC_SITE_DIST          => 19,
    JSRC_ALN_FAILURE_FLAG   => 20,
    JSRC_JXN_FAILURE_FLAG   => 21,
    JSRC_SEQS               => 22,
    JSRC_QUALS              => 23,
    JSRC_CIGARS             => 24,
};

# open file handles
open my $jxnH, "-|", "zcat $SV_UNIQUE_JUNCTIONS_FILE" or die "could not open $SV_UNIQUE_JUNCTIONS_FILE : $!\n";
open my $srcH, "-|", "zcat $SV_JUNCTION_SOURCES_FILE" or die "could not open $SV_JUNCTION_SOURCES_FILE : $!\n";

# process junctions one key at a time
my $prevSrcLine = <$srcH>;
while (my $jxn = <$jxnH>){
    chomp $jxn;
    my @jxn = split("\t", $jxn);
    my $jxnKeyWithoutJxnProp = join("\t", @jxn[OJXN_CHROM_INDEX1_1..OJXN_JXN_TYPE]);
    my $jxnKeyWIthJxnProp    = join("\t", $jxnKeyWithoutJxnProp, @jxn[UJXN_ALN_OFFSET..UJXN_JXN_BASES]);
    while($prevSrcLine and $prevSrcLine =~ m/^$jxnKeyWithoutJxnProp/){
        chomp $prevSrcLine;
        my @src = split("\t", $prevSrcLine);

        # values taken from the best, most frequent junction signature
        $prevSrcLine =~ m/^$jxnKeyWIthJxnProp/ and 
            @jxn[UJXN_JSRC_MAPQ..UJXN_JSRC_JXN_FAILURE_FLAG] = @src[JSRC_MAPQ..JSRC_JXN_FAILURE_FLAG]; 

        # properties collected across all junction keys regardles of junction properties for PCR duplciation/ONT duplex purging
        foreach my $pair(
            [UJXN_JSRC_QNAMES,           JSRC_QNAMES], # all fields made comma-safe upstream
            [UJXN_JSRC_INSERT_SIZES,     JSRC_INSERT_SIZES],
            [UJXN_JSRC_IS_ALLOWED_SIZES, JSRC_IS_ALLOWED_SIZES],
            [UJXN_JSRC_MIN_STEM_LENGTHS, JSRC_MIN_STEM_LENGTHS],
            [UJXN_JSRC_PASSED_STEMS,     JSRC_PASSED_STEMS],
            [UJXN_JSRC_SEQS,             JSRC_SEQS],
            [UJXN_JSRC_QUALS,            JSRC_QUALS],
            [UJXN_JSRC_CIGARS,           JSRC_CIGARS],
            [UJXN_JSRC_ORIENTATIONS,     JSRC_ORIENTATIONS],
            [UJXN_JSRC_CHANNELS,         JSRC_CHANNELS],
            [UJXN_JSRC_OUTER_ENDPOINTS,  JSRC_OUTER_ENDPOINTS]
        ) {
            my $jxnI = $$pair[0];
            my $srcI = $$pair[1];
            $jxn[$jxnI] = 
                defined $jxn[$jxnI] ? 
                join(",", $jxn[$jxnI], $src[$srcI]) : 
                $src[$srcI];
        }
        $prevSrcLine = <$srcH>;
    }
    print join("\t", @jxn), "\n";
}

# close file handles
close $jxnH;
close $srcH;

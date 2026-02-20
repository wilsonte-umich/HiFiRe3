use strict;
use warnings;

# action: create a composite genome metadata file for ues in the app Track Browser

my $GENOME1 = $ENV{GENOME1};
my $GENOME2 = $ENV{GENOME2};
my $GENOME  = $GENOME1 . "_" . $GENOME2;

my $templateFile = $ENV{ACTION_DIR} . "/app_metadata_template.yml";
my $metadataFile = $ENV{GENOME_DIR} . "/$GENOME.yml";

# read template into string
open my $fh, '<', $templateFile or die "Could not open file '$templateFile' $!";
my $template = do { local $/; <$fh> };
close $fh;

# set the gene prediction track types based on the genome names (hs1 requires special handling)
my $genePred1 = $GENOME1 eq "hs1" ? "bigGenePred" : "genePred";
my $genePred2 = $GENOME2 eq "hs1" ? "bigGenePred" : "genePred";
my $track1    = $GENOME1 eq "hs1" ? "hub_3671779_ncbiRefSeqCurated" : "ncbiRefSeqCurated";
my $track2    = $GENOME2 eq "hs1" ? "hub_3671779_ncbiRefSeqCurated" : "ncbiRefSeqCurated";

# replace placeholders with actual values
$template =~ s/__GENOME_1__/$GENOME1/g;
$template =~ s/__GENOME_2__/$GENOME2/g;
$template =~ s/__GENE_PRED_1__/$genePred1/g;
$template =~ s/__GENE_PRED_2__/$genePred2/g;
$template =~ s/__TRACK_1__/$track1/g;
$template =~ s/__TRACK_2__/$track2/g;

# write the modified template to the metadata file
open my $out, '>', $metadataFile or die "Could not open file '$metadataFile' $!";
print $out $template;
close $out;

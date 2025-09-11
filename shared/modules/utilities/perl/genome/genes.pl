use strict;
use warnings;

#----------------------------------------------------------
# manipulations related to gene correspondence of breakpoints
#----------------------------------------------------------

# output variables
our ($nGenes, $sumGeneLens) = (0, 0);
our (%genes, @geneInfo);

# constants
use constant {
    NULL_GENE => ",",
    NULL_GENE_DISTANCE => ",",
};

# operating parameters
use vars qw($GENES_BED $GENE_SCALAR);
defined $$GENE_SCALAR or $$GENE_SCALAR  = 1;

# load the genes
sub loadGeneRegions {
    my ($quiet) = @_;
    ($GENES_BED and $GENES_BED ne "null" and $GENES_BED ne "NA") or return;
    @geneInfo = ();
    $nGenes = loadGeneRegions_(\$sumGeneLens);
    if(!$quiet){
        printCount(commify($nGenes),      'nGenes',      'genes');
        printCount(commify($sumGeneLens), 'sumGeneLens', 'total bp covered by genes');
    }
}
sub loadGeneRegions_ {
    my ($sumLens) = @_;
    open my $inH, "-|", "zcat $GENES_BED" or die "could not open $GENES_BED: $!\n";
    my $geneI1 = 0;
    while(my $line = <$inH>){
        $geneI1++; # 1-referenced
        chomp $line;
        $line =~ s/\r//g;
        my ($chr, $start, $end, $name) = split("\t", $line);
        for(my $pos =  int($start / $GENE_SCALAR);
               $pos <= int($end   / $GENE_SCALAR);
               $pos++){
            push @{$genes{$chr}{$pos}}, $geneI1; # lookup for the gene associated with a coordinate
        }
        my $len = $end - $start;
        $$sumLens += $len;
        push @geneInfo, [$start + int($len / 2), $name];
    }
    close $inH;
    return $geneI1;
}

# return information on a single genome position relative to genes
sub getPosGene {
    my ($chr, $pos) = @_;
    my $cg = $genes{$chr} or return (NULL_GENE, NULL_GENE_DISTANCE);
    my $cgp = $$cg{int($pos / $GENE_SCALAR)} || 0;
    if($cgp){
        my $genes = NULL_GENE;
        my $dists = NULL_GENE_DISTANCE;
        foreach my $geneI1 (@$cgp){
            my ($geneCenter, $name) = @{$geneInfo[$geneI1 - 1]};
            $genes .= "$name,";
            $dists .= $pos - $geneCenter . ",";
        }
        ($genes, $dists);
    } else {
        (NULL_GENE, NULL_GENE_DISTANCE);
    }
}

1;

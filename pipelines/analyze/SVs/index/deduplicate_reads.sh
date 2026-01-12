# action:
#   sort and count bases in a library
#   record only one value per unique path alignment to deduplicate counts
# input:
#   alignment_sizes files over all chromosomes on STDIN (since chrom-level files have re-ordered translocation rows)
# output:
#   aggregated alignments on STDOUT

# orderedNode1,orderedNode2,channel,(reordered)alnI0,nRefBases,nReadBases
$SORT -k1,1n -k2,2n -k3,3n -k4,4n |
bedtools groupby -g 1,2,3,4 -c 5,6 -o max,max # deduplication occurs here

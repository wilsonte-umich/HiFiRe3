# action:
#   sort and count bases in a library
#   record only one value per unique path alignment to deduplicate counts
# input:
#   alignment_sizes files over all chromosomes on STDIN (since chrom-level files have re-ordered translocation rows)
# output:
#   aggregated alignments on STDOUT

# chromIndex1_1,siteIndex1_1,chromIndex1_2,siteIndex1_2,channel,(reordered)alnI0,sumRefBases,sumReadBases # ,end5OnTarget
$SORT -k1,1n -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n |
bedtools groupby -g 1,2,3,4,5,6 -c 7,8 -o max,max | # deduplication occurs here
bedtools groupby -g 1,3,6 -c 7,8,8 -o sum,sum,count 
# chromIndex1_1,chromIndex1_2,(reordered)alnI0,sumRefBases,sumReadBases,nAlns # ,end5OnTarget

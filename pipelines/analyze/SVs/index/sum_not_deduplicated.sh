# action:
#   sort and count bases in a library
#   sum all values without deduplication
# input:
#   alignment_sizes files over all chromosomes on STDIN (since chrom-level files have re-ordered translocation rows)
# output:
#   aggregated alignments on STDOUT

# chromIndex1_1,chromIndex1_2,(reordered)alnI0,sumRefBases,sumReadBases # ,end5OnTarget
$SORT -k1,1n -k2,2n -k3,3n |
bedtools groupby -g 1,2,3 -c 4,5,5 -o sum,sum,count 
# chromIndex1_1,chromIndex1_2,(reordered)alnI0,sumRefBases,sumReadBases,nAlns # ,end5OnTarget

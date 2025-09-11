# actions:
#   report some useful summary information to answer questions like:
#       what fraction of the genome was on target, in genes, both, etc.?
# input:
#   $TARGETS_BED
#   $GENES_BED
#   $GENOME_FASTA.fai
# output:
#   BED file of genome split by target and gene region states
#   count of bases and segments in each state combination

export GROUP_BY=TRUE
export SCORE_TYPE=COUNT
export NAME_JOIN_CHAR=,
export INCLUDE_GAPS=TRUE
export CHROM_FILE=${GENOME_FASTA}.fai
export USE_NONCANONICAL=TRUE # to ensure that concatenated chromosomes are retained

BED_UTIL_SPLIT=${MODULES_DIR}/utilities/perl/genome/bedutil/split.pl

echo "tallying genome regions by target and gene states"

cat <(
    if [ -f $TARGETS_BED ]; then
        cat $TARGETS_BED | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "TARGET_"$4, 0, "+"}'
    fi
) <(
    zcat $GENES_BED | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "GENE_"$4, 0, "+"}'
) | 
perl ${BED_UTIL_SPLIT} |
sort -k1,1 -k2,2n -k3,3n --parallel=$N_CPU --buffer-size=4G |
grep -v -P "\w+_\w+_\w+" | 
perl ${ACTION_DIR}/tally/tally_genome.pl |
tee ${GENOME_TALLY_FILE}
checkPipe

echo "done"

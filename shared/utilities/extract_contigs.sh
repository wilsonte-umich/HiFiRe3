# action:
#   create a FASTA file of the non-gapped, non-blacklist contigs from a reference genome
#   to be used as a "haplotype" along with the same reference to build a pseudograph assembly
#   this is a useful approach for inbred mouse strains, etc.
#   the pseudograph is used in genotype analysis in place of a true diploid graph assembly
# expects:
#   $GENOME_DIR, as populated by `prepare genome`
#   bedtools
# output:
#   FASTA stream on stdout
# usage:
#   module load bedtools # as needed
#   export GENOME_DIR=/path/to/prepared/genome
#   bash extract_contigs.sh > /path/to/haplotype.fa

# collect the genome
if [ "$GENOME_DIR" = "" ]; then
    echo "missing environment variable: GENOME_DIR"
    exit 1
fi
GENOME=`basename $GENOME_DIR`
GENOME_FASTA=${GENOME_DIR}/${GENOME}.fa
GENOME_INDEX=${GENOME_FASTA}.fai

# get the exclusion files
METADATA_DIR=${GENOME_DIR}/metadata
GAPS_FILE=${METADATA_DIR}/${GENOME}.gaps.txt
EXCLUSIONS_FILE=${METADATA_DIR}/${GENOME}.exclusions.bed

# prepare the contig regions as a BED stream
cat <(
    if [ -f ${GAPS_FILE} ]; then
        cut -f 2-4 ${GAPS_FILE}
    else 
        cat /dev/null
    fi
) <(
    if [ -f ${EXCLUSIONS_FILE} ]; then
        cut -f 1-3 ${EXCLUSIONS_FILE}
    else 
        cat /dev/null
    fi
) |
sort -k1,1 -k2,2n -k3,3n | 
bedtools complement -i - -g ${GENOME_INDEX} | # thus, BED now includes only the non-excluded spans
grep -v '_' | # on canonical chromosomes only

# extract the corresponding genome sub-sequences to FASTA
bedtools getfasta -fi ${GENOME_FASTA} -bed - | 
fold -w 60

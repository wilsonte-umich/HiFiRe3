# action:
#     set environment variables for genome file access
#     index $GENOME.fa as needed
#     collect chromosome metadata
# expects:
#     $GENOME_DIR, will default to $MDI_DIR/resources/genomes/$GENOME
#     $GENOME
#     indexed fasta file $GENOME.fa or genome.fa in $GENOME_DIR
# usage:
#     source $MODULES_DIR/genome/set_genome_vars.sh

# file prefixes for genome-specific output files
export DATA_GENOME_PREFIX=${DATA_FILE_PREFIX}.${GENOME}
export PLOT_GENOME_PREFIX=${PLOT_PREFIX}.${GENOME}

# set and check the genome directory
if [[ "$GENOME_DIR" == "" || "$GENOME_DIR" == "null" || "$GENOME_DIR" == "NA" ]]; then
    export GENOME_DIR=${MDI_DIR}/resources/genomes/${GENOME}
    mkdir -p ${GENOME_DIR}
fi
if [ ! -d $GENOME_DIR ]; then
    echo
    echo "--genome-dir not found:"
    echo $GENOME_DIR
    echo
    exit 1
fi
export GENOME_PREFIX=${GENOME_DIR}/${GENOME}
export GENOME_FASTA=${GENOME_PREFIX}.fa

# metadata directories and files
export GENOME_METADATA_DIR=${GENOME_DIR}/metadata
export GENOME_METADATA_PREFIX=${GENOME_METADATA_DIR}/${GENOME}
export GENOME_GAPS_FILE=${GENOME_METADATA_PREFIX}.gaps.txt
export GENOME_EXCLUSIONS_BED=${GENOME_METADATA_PREFIX}.exclusions.bed

# annotations
export GENOME_ANNOTATIONS_DIR=${GENOME_DIR}/annotations
export ANNOTATION_GTF=${GENOME_ANNOTATIONS_DIR}/${GENOME}.ncbiRefSeq.gtf.gz
export GENES_BED=${GENOME_ANNOTATIONS_DIR}/${GENOME}.ncbiRefSeq.genes.bed.gz

# initialize genome directory tree
if [ "$HIFIRE3_PREPARING_GENOME" == "TRUE" ]; then
    mkdir -p ${GENOME_METADATA_DIR}
    mkdir -p ${GENOME_ANNOTATIONS_DIR}

# check and collect genome metadata
else 
    # fasta file and index
    if [ ! -f $GENOME_FASTA ]; then
        export GENOME_FASTA=${GENOME_DIR}/genome.fa
    fi
    if [ ! -f $GENOME_FASTA ]; then
        echo "missing genome fasta file in ${GENOME_DIR}"
        echo "expected either file ${GENOME}.fa or genome.fa"
        exit 1
    fi
    if [ ! -f $GENOME_FASTA.fai ]; then
        echo "indexing genome fasta file"
        samtools index ${GENOME_FASTA}
    fi

    # get genome size
    # get the list of all placed chromosome sequences, including chrX, chrY, chrM, and chrEBV if present
    # if instructed, use all chromosomes regardlesss of name 
    if [[ "$USE_ALL_CHROMS" != "" &&  "$USE_ALL_CHROMS" != "0" ]]; then
        export GENOME_SIZE=`awk '{s+=$2}END{print s}' $GENOME_FASTA.fai`
        export GENOME_CHROMS=`cat $GENOME_FASTA.fai | cut -f1`

    # if a composite spike-in genome, filter against non-canonical chroms with two _ characters in name
    # since we expect canonical chrom names in format <source chrom name>_<source genome name>, e.g. chr1_hg38
    elif [[ "$IS_COMPOSITE_GENOME" != "" &&  "$IS_COMPOSITE_GENOME" != "0" ]]; then
        export GENOME_SIZE=`awk '$1!~/_.*_/ && $1!="chrEBV"{s+=$2}END{print s}' $GENOME_FASTA.fai`
        export GENOME_CHROMS=`cat $GENOME_FASTA.fai | cut -f1 | grep -v '_.*_'`

    # otherwise, filter against non-canonical chroms with any _ character in name
    else 
        export GENOME_SIZE=`awk '$1!~/_/ && $1!="chrEBV"{s+=$2}END{print s}' $GENOME_FASTA.fai`
        export GENOME_CHROMS=`cat $GENOME_FASTA.fai | cut -f1 | grep -v _`
    fi
fi

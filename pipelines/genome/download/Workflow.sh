#!/bin/bash

show_alert_message () {
    echo
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo -e "$1"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo
}

# check that the requested genome is supported
export SUPPORTED_GENOMES="hs1|CHM13|hg38|GRCh38|mm39|GRCm39|dm6"
if [[ "|${SUPPORTED_GENOMES}|" != *"|${GENOME}|"* ]]; then
    show_alert_message \
"genome not supported for download: $GENOME\n"\
"supported genomes names are: $SUPPORTED_GENOMES"
    exit 1
fi

# coerce GRC to UCSC genome names
# using UCSC names helps support the MDI genome browser
#   which dynamically reads genome metadata from the UCSC API
export UCSC_GENOMES="hs1|hg38|mm39|dm6"
if [[ "|${UCSC_GENOMES}|" != *"|${GENOME}|"* ]]; then
    GRC_GENOME=${GENOME}
    if [ "$GENOME" = "CHM13" ]; then
        export GENOME=hs1
    elif [ "$GENOME" = "GRCh38" ]; then
        export GENOME=hg38
    elif [ "$GENOME" = "GRCm39" ]; then
        export GENOME=mm39
    fi
    show_alert_message \
"GRC genome name ${GRC_GENOME} was coerced to the UCSC genome name ${GENOME}\n"\
"usage is optimized for UCSC genome names as they directly support the genome browser"
fi

# set derivative environment variables and file paths
# don't use set_genome_vars.sh, we don't want its cascading actions here
#------------------------------------------------------------------------
export GENOME_DIR=${GENOMES_DIR}/${GENOME}
#------------------------------------------------------------------------
export GENOME_METADATA_DIR=${GENOME_DIR}/metadata
mkdir -p ${GENOME_METADATA_DIR}
export GENOME_METADATA_PREFIX=${GENOME_METADATA_DIR}/${GENOME}
export GENOME_GAPS_FILE=${GENOME_METADATA_PREFIX}.gaps.txt
export GENOME_EXCLUSIONS_BED=${GENOME_METADATA_PREFIX}.exclusions.bed
#------------------------------------------------------------------------
export GENOME_ANNOTATIONS_DIR=${GENOME_DIR}/annotations
export ANNOTATION_GTF=${GENOME_ANNOTATIONS_DIR}/${GENOME}.ncbiRefSeq.gtf.gz
export GENES_BED=${GENOME_ANNOTATIONS_DIR}/${GENOME}.ncbiRefSeq.genes.bed.gz
mkdir -p ${GENOME_ANNOTATIONS_DIR}
#------------------------------------------------------------------------
export GENOME_PREFIX=${GENOME_DIR}/${GENOME}
export GENOME_FASTA=${GENOME_PREFIX}.fa

# download the genome fasta and supporting metadata
runWorkflowStep 1 download $ACTION_DIR/download.sh

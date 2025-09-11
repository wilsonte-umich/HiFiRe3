#!/bin/bash

# set primary genome paths
export GENOME1_DIR=${GENOMES_DIR}/${GENOME1}
#------------------------------------------------------------------------
export GENOME1_METADATA_DIR=${GENOME1_DIR}/metadata
export GENOME1_METADATA_PREFIX=${GENOME1_METADATA_DIR}/${GENOME1}
export GENOME1_GAPS_FILE=${GENOME1_METADATA_PREFIX}.gaps.txt
export GENOME1_EXCLUSIONS_BED=${GENOME1_METADATA_PREFIX}.exclusions.bed
#------------------------------------------------------------------------
export GENOME1_ANNOTATIONS_DIR=${GENOME1_DIR}/annotations
export GENES1_BED=${GENOME1_ANNOTATIONS_DIR}/${GENOME1}.ncbiRefSeq.genes.bed.gz
#------------------------------------------------------------------------
export GENOME1_PREFIX=${GENOME1_DIR}/${GENOME1}
export GENOME1_FASTA=${GENOME1_PREFIX}.fa

# set secondary genome paths
export GENOME2_DIR=${GENOMES_DIR}/${GENOME2}
#------------------------------------------------------------------------
export GENOME2_METADATA_DIR=${GENOME2_DIR}/metadata
export GENOME2_METADATA_PREFIX=${GENOME2_METADATA_DIR}/${GENOME2}
export GENOME2_GAPS_FILE=${GENOME2_METADATA_PREFIX}.gaps.txt
export GENOME2_EXCLUSIONS_BED=${GENOME2_METADATA_PREFIX}.exclusions.bed
#------------------------------------------------------------------------
export GENOME2_ANNOTATIONS_DIR=${GENOME2_DIR}/annotations
export GENES2_BED=${GENOME2_ANNOTATIONS_DIR}/${GENOME2}.ncbiRefSeq.genes.bed.gz
#------------------------------------------------------------------------
export GENOME2_PREFIX=${GENOME2_DIR}/${GENOME2}
export GENOME2_FASTA=${GENOME2_PREFIX}.fa

# set output concatenated genome paths
#------------------------------------------------------------------------
export GENOME=${GENOME1}_${GENOME2}
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
runWorkflowStep 1 concatenate $ACTION_DIR/concatenate.sh

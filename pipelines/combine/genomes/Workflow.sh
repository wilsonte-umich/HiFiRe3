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

# set derivative environment variables and file paths
export GENOME=${GENOME1}_${GENOME2}
export GENOME_DIR=${GENOMES_DIR}/${GENOME}
mkdir -p ${GENOME_DIR}
export HIFIRE3_PREPARING_GENOME=TRUE
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/REs/set_re_vars.sh

# download the genome fasta and supporting metadata
runWorkflowStep 1 concatenate concatenate.sh

# update derivative environment variables with genome metadata
export HIFIRE3_PREPARING_GENOME=FALSE
export IS_COMPOSITE_GENOME=TRUE
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/REs/set_re_vars.sh

# collect all genome RE sites for a set of restriction enzymes
runWorkflowStep 2 digest ${MODULES_DIR}/genome/digest.sh

# analyze the restriction fragment length distributions for plotting in app
runWorkflowStep 3 assemble ${MODULES_DIR}/genome/assemble.sh

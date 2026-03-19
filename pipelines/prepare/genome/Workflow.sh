#!/bin/bash

# parse requested genome
source ${MODULES_DIR}/genome/parse_genome.sh

# set derivative environment variables and file paths
export GENOME_DIR=${GENOMES_DIR}/${GENOME}
mkdir -p ${GENOME_DIR}
export HIFIRE3_PREPARING_GENOME=TRUE
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/REs/set_re_vars.sh

# download the genome fasta and supporting metadata
runWorkflowStep 1 download ${MODULES_DIR}/genome/download.sh

# update derivative environment variables with genome metadata
export HIFIRE3_PREPARING_GENOME=FALSE
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/REs/set_re_vars.sh

# collect all genome RE sites for a set of restriction enzymes
runWorkflowStep 2 digest ${MODULES_DIR}/genome/digest.sh

# analyze the restriction fragment length distributions for plotting in app
runWorkflowStep 3 assemble ${MODULES_DIR}/genome/assemble.sh

#!/bin/bash

# check the genome installation and prepare it if necessary
runWorkflowStep 1 check_genome $MODULES_DIR/genome/check_genome.sh

# # set derivative environment variables and file paths
source $MODULES_DIR/rust/set_rust_vars.sh
source $MODULES_DIR/genome/set_genome_vars.sh
UBAM_DIR=${TASK_DIR}/ubam
mkdir -p ${UBAM_DIR}
export READ_BAM_FILE=${UBAM_DIR}/${DATA_NAME}.reads.unaligned.bam
export PACBIO_BASECALL_KINETICS=${UBAM_DIR}/${DATA_NAME}.kinetics.unaligned.txt.gz

# convert PacBio strand-specific BAM files to read consensus BAM with kinetics
runWorkflowStep 2 basecall $ACTION_DIR/basecall.sh

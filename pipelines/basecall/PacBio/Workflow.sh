#!/bin/bash

# # set derivative environment variables and file paths
source $MODULES_DIR/rust/set_rust_vars.sh
source $MODULES_DIR/genome/set_genome_vars.sh
UBAM_DIR=${TASK_DIR}/ubam
mkdir -p ${UBAM_DIR}
export READ_BAM_FILE=${UBAM_DIR}/${DATA_NAME}.reads.unaligned.bam
export PACBIO_BASECALL_KINETICS=${UBAM_DIR}/${DATA_NAME}.kinetics.unaligned.txt.gz

# convert PacBio strand-specific BAM files to read consensus BAM with kinetics
runWorkflowStep 1 basecall $ACTION_DIR/basecall.sh

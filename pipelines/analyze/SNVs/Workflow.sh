#!/bin/bash

# check the genome installation and prepare it if necessary
runWorkflowStep 1 check_genome $MODULES_DIR/genome/check_genome.sh

# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/workflow.sh
export NAME_BAM_FILES=${NAME_BAM_FILE}

# index RE fragments for unique SNV/indel combinations
# create a pileup of alignments relative to genome coordinates
# call (sub)clonal SNVs and indels
runWorkflowStep 2 analyze_SNVs $MODULES_DIR/analyze/analyze_SNVs.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM

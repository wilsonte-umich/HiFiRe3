#!/bin/bash

# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/workflow.sh
export NAME_BAM_FILES=${NAME_BAM_FILE}

# index RE fragments for unique SNV/indel combinations
# create a pileup of alignments relative to genome coordinates
# call (sub)clonal SNVs and indels
runWorkflowStep 1 analyze_SNVs $MODULES_DIR/analyze/analyze_SNVs.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM

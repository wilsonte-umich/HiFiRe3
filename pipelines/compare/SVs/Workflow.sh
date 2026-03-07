#!/bin/bash

# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/workflow.sh

# group and tally junctions across multiple samples
runWorkflowStep 1 compare_SNVs $MODULES_DIR/analyze/analyze_SVs.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM

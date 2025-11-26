#!/bin/bash

# set derivative environment variables and file paths
export VARIANT_TYPE=SVs
export VARIANT_FILE_SUFFIX=final_junctions_1.txt.bgz
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/workflow.sh

# group and tally junctions across multiple samples
runWorkflowStep 1 compare compare.sh

# reduce the size of the final junctions files for loading in the app
runWorkflowStep 2 package ${MODULES_DIR}/analyze/SVs/package/package.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM

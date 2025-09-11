#!/bin/bash

# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/Workflow_1.sh
mkdir -p $PLOTS_DIR

# match alignments to RE sites and fragments
# apply various SV quality filters and error correction mechanisms
# write SITE_SAM for subsequent fragment and variant aggregation
runWorkflowStep 1 match_sites match_sites/match_sites.sh

# plot the insert size distribution and calculate insert size representation in library
runWorkflowStep 2 insert_sizes insert_sizes/insert_sizes.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM

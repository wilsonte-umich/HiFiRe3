#!/bin/bash

# check the genome installation and prepare it if necessary
runWorkflowStep 1 check_genome $MODULES_DIR/genome/check_genome.sh

# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/workflow.sh
export NAME_BAM_FILES=${NAME_BAM_FILE}

# index RE fragments for unique SV paths
# perform fuzzy matching of junction nodes to each other to aggregate inexact junction matches
# purge ONT duplexes when junctions appeared on different strands in the same channel
# record number of supporting reads per junction and whether any/all were inside/outside the allowed insert size range
# create a second reverse-sorted junction file for indexed retrieval of both junction nodes for independent local plotting
runWorkflowStep 2 analyze_SVs $MODULES_DIR/analyze/analyze_SVs.sh

# # report some useful summary information for junctions
# runWorkflowStep 3 tally_junctions tally/tally_junctions.sh

# # report some useful summary information target/genic genome regions
# runWorkflowStep 4 tally_genome tally/tally_genome.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM

#!/bin/bash

# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/Workflow.sh
source $MODULES_DIR/align/set_read_file_vars.sh
mkdir -p $PLOTS_DIR

#-------------------------------------------------------------------------------
# read alignment
#-------------------------------------------------------------------------------

# set working directory to READ_FILE_DIR to avoid too-long argument list with multiple read files
cd ${READ_FILE_DIR}

# align read sequences to reference genome
runWorkflowStep 1 align align/align.sh

# reset working directory
cd ${TASK_DIR}

#-------------------------------------------------------------------------------
# RE site localization
#-------------------------------------------------------------------------------

# make a tally of the informative endpoints of all reads
# many/most should match in silico RE sites below
runWorkflowStep 2 extract_endpoints locate/extract_endpoints.sh

# match endpoints to each other and to in silico site positions
# creates a table of filtering sites for tolerance matching, etc.
runWorkflowStep 3 tabulate_endpoints locate/tabulate_endpoints.sh

# TODO: consider letting match_sites.sh create the index files in temporary directories
#       rather than writing them to the permanent output directory

# create binary lookup files to speed matching of sample endpoints to filtering sites
runWorkflowStep 4 create_index locate/create_index.sh

#-------------------------------------------------------------------------------
# read alignment parsing, including matching to RE sites
#-------------------------------------------------------------------------------

# match alignments to RE sites and fragments
# apply various alignment and SV quality filters and error correction mechanisms
# write SITE_SAM for subsequent fragment and variant aggregation
runWorkflowStep 5 apply_filters apply_filters/apply_filters.sh

#-------------------------------------------------------------------------------
# insert size analysis
#-------------------------------------------------------------------------------

# create plots of insert size distributions before and after filtering and RE site projection
# establish automated thresholds for allowed insert sizes working with the hint from --min-selected-size
runWorkflowStep 6 insert_sizes insert_sizes/analyze_insert_sizes.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM

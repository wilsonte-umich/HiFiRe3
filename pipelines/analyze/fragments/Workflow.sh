#!/bin/bash

# check the genome installation and prepare it if necessary
runWorkflowStep 1 check_genome $MODULES_DIR/genome/check_genome.sh

# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/workflow.sh
source $MODULES_DIR/align/set_read_file_vars.sh
mkdir -p $PLOTS_DIR
export ALIGN_DIR=${ACTION_DIR}/align

#-------------------------------------------------------------------------------
# read alignment, initial fragment analysis, and RE endpoint extraction
#-------------------------------------------------------------------------------

# set working directory to READ_FILE_DIR to avoid too-long argument list with multiple read files
cd ${READ_FILE_DIR}

# align read sequences to reference genome, or just analyze Ultima aligned reads
runWorkflowStep 2 align align/align.sh

# reset working directory
cd ${TASK_DIR}

#-------------------------------------------------------------------------------
# RE site localization
#-------------------------------------------------------------------------------

# match endpoints to each other and to in silico site positions
# create a table of filtering sites for tolerance matching, etc.
runWorkflowStep 3 tabulate_endpoints locate/tabulate_endpoints.sh

#-------------------------------------------------------------------------------
# read alignment parsing, including matching to RE sites
#-------------------------------------------------------------------------------

# match alignments to RE sites and fragments to characterize inserts
# apply various alignment and SV quality filters and error correction mechanisms
runWorkflowStep 4 analyze_inserts analyze_inserts.sh

#-------------------------------------------------------------------------------
# insert size analysis
#-------------------------------------------------------------------------------

# create plots of insert size distributions before and after filtering and RE site projection
# establish automated thresholds for allowed insert sizes working with the hint from --min-selected-size
runWorkflowStep 5 insert_sizes insert_sizes/analyze_insert_sizes.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM

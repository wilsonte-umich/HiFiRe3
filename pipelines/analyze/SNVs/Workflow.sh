#!/bin/bash

# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/Workflow_1.sh
export TMP_PILEUP_DIR=$TMP_DIR_WRK_SMALL

# index RE fragments for unique SNV/indel combinations
runWorkflowStep 1 index_fragments index/index_fragments.sh

# create a pileup of alignments relative to genome coordinates
# call (sub)clonal SNVs and indels
# create a pileup of alignments relative to read positions
runWorkflowStep 2 pileup pileup/pileup.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM
rm -f $SNV_GENOME_PILEUP_PREFIX.*.bed.gz
rm -f $SNV_SUMMARY_TABLE_PREFIX.*.bed.gz
rm -f $SNV_READ_PILEUP_PREFIX.*.bed.gz

#!/bin/bash

# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/Workflow.sh

# copy GENOME_FASTA to TMP_DIR_WRK_SHM for fast access
echo "copying $GENOME fasta file to shared memory"
cp $GENOME_FASTA     $TMP_DIR_WRK_SHM
cp $GENOME_FASTA.fai $TMP_DIR_WRK_SHM
export GENOME_FASTA_SHM=$TMP_DIR_WRK_SHM/$(basename $GENOME_FASTA) 

# index RE fragments for unique SV paths
runWorkflowStep 1 index_fragments index/index_fragments.sh

# merge event metadata from junction_sources into junctions
# perform fuzzy matching of junction nodes to each other to aggregate inexact junction matches
# purge ONT duplexes when junctions appeared on different strands in the same channel
# create a second reverse-sorted junction file for indexed retrieval of both junction nodes for independent local plotting
runWorkflowStep 2 aggregate_SVs aggregate/aggregate_SVs.sh

# report some useful summary information for junctions
runWorkflowStep 3 tally_junctions tally/tally_junctions.sh

# report some useful summary information target/genic genome regions
runWorkflowStep 4 tally_genome tally/tally_genome.sh

# report some useful summary information target/genic genome regions
runWorkflowStep 5 package package/package.sh

# clean up
rm -fr $TMP_DIR_WRK_SMALL
rm -fr $TMP_DIR_WRK_SHM

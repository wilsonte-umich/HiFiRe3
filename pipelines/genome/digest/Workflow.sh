#!/bin/bash

# action 'genome digest' must be run once by all users
# to initialize RE enzyme cleavage sites expected for $GENOME

# set derivative environment variables and file paths
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/REs/set_re_vars.sh
mkdir -p $GENOME_REMAPS_DIR

# set output paths
export RE_SUMMARIES_RDS=${GENOME_PREFIX}.RE_site_summaries.rds
export RE_DISTRIBUTIONS_RDS=${GENOME_PREFIX}.RE_distributions.rds
export RE_FRAGMENTS_RDS=${GENOME_PREFIX}.RE_fragments.rds

# collect all genome RE sites for a set of restriction enzymes
runWorkflowStep 1 digest digest.sh

# analyze the restriction fragment length distributions for plotting in app
runWorkflowStep 2 assemble assemble.sh

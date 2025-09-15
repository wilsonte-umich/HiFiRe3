#!/bin/bash


# TODO:
# for HiFiRe3, we will only score SNVs/indels for:
#     PacBio reads (not ONT), an end-to-end platform
#     RE fragments with sufficient coverage depth
# therefore, need to:
#     filter alignments to non-SV molecules
#     sort alignments by sitePos1,sitePos2
#     extract groups of alignments with the same sitePos1,sitePos2 (both strands)
#     create a pileup of alignments for each RE fragment group
# need CS_TAG to be carried through to SITE_SAM

# current code uses the va tag in a way that preserved the order on the read
# motivated by a decay in accuracy of short reads over the read
# this decay isn't systematically expected for HiFi reads, so stop doing this most likely
# probably just need to track variants relative to the top genome strand, as cs tag does
# HiFi reads are also largely unstranded EXCEPT for heterouplex
# so will want to carry the strand along 


# set derivative environment variables and file paths
export PIPELINE_SHARED_DIR=${PIPELINE_DIR}/shared
source ${PIPELINE_SHARED_DIR}/Workflow.sh
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

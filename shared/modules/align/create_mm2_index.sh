#!/bin/bash

# actions:
#     as needed, create and save an appropriate minimap2 genome alignment index
# input:
#     ${GENOME_FASTA_WRK} set to the appropriate path by caller
#     ${ALIGNMENT_MODE_WRK}
# output:
#     ${MINIMAP2_INDEX_WRK}, variable set by this script based on above

# set the path to the minimap2 index, name is a suffixed version of ${GENOME_FASTA_WRK}
export MINIMAP2_INDEX_WRK=${GENOME_FASTA_WRK}.${ALIGNMENT_MODE_WRK}.mmi

# as needed, create and save an appropriate minimap2 genome alignment index
# index is stored for future use also
# use a lock to prevent concurrent jobs from trying to create the same index
MINIMAP2_INDEX_LOCK=${MINIMAP2_INDEX_WRK}.lock
while [ -f ${MINIMAP2_INDEX_LOCK} ]; do # prevent two jobs from creating the same index at the same time
    echo "waiting for minimap2 index to be created"
    sleep 60
done
if [ ! -f "${MINIMAP2_INDEX_WRK}" ]; then
    touch ${MINIMAP2_INDEX_LOCK}
    echo
    echo "creating minimap2 index for:"
    echo "    fasta:  "`basename ${GENOME_FASTA_WRK}`
    echo "    mode:   "${ALIGNMENT_MODE_WRK}
    minimap2 -x ${ALIGNMENT_MODE_WRK} -t 3 -d ${MINIMAP2_INDEX_WRK} ${GENOME_FASTA_WRK}
    checkPipe
    rm ${MINIMAP2_INDEX_LOCK}
    echo "minimap2 index done"
    if [ "$CREATE_MM2_INDEX_MESSAGE" != "" ]; then
        echo
        echo $CREATE_MM2_INDEX_MESSAGE
    fi
fi

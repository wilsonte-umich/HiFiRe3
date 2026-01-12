#!/bin/bash

# working variables
export SORT_COMMAND="sort --parallel=$N_CPU --buffer-size=4G --temporary-directory $TMP_DIR_WRK_SMALL"
BGZIP="bgzip --threads $N_CPU --force --output"
TABIX="tabix --threads $N_CPU"

# sort the junctions by chromIndex1_1,strandIndex0_1,refPos1_1 in preparation for breakpoint node1 fuzzy matching
zcat $SAMPLE_FILES | 
$SORT_COMMAND -k1,1n -k3,3n -k2,2n |
perl ${ACTION_DIR}/compare.pl | 

# sort and index the junctions by node1
$SORT_COMMAND -k1,1n -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n|
$BGZIP $SV_FINAL_JUNCTIONS_FILE_1
checkPipe
$TABIX --sequence 1 --begin 2 --end 2 $SV_FINAL_JUNCTIONS_FILE_1
checkPipe

# resort and index the junctions file by node2
zcat $SV_FINAL_JUNCTIONS_FILE_1 | 
$SORT_COMMAND -k4,4n -k5,5n -k6,6n -k1,1n -k2,2n -k3,3n | 
$BGZIP $SV_FINAL_JUNCTIONS_FILE_2
checkPipe
$TABIX --sequence 4 --begin 5 --end 5 $SV_FINAL_JUNCTIONS_FILE_2
checkPipe

echo
echo "done"

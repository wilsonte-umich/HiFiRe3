#!/bin/bash

# pack index files for fast retrieval of the RE site positions closest to endpoints

BGZIP="bgzip --threads $N_CPU --force --output"
TABIX="tabix --threads $N_CPU"

echo "creating site lookup files"
perl ${ACTION_DIR}/locate/create_index.pl | 
$BGZIP $FILTERING_SITES_BGZ
checkPipe

echo "creating tabix index"
$TABIX --sequence 1 --begin 2 --end 2 $FILTERING_SITES_BGZ
checkPipe

echo "done"

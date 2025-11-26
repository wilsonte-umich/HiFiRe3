#!/bin/bash

# pack index files for fast retrieval of the RE site positions closest to endpoints
BGZIP="bgzip --threads $N_CPU --force --output"
TABIX="tabix --threads $N_CPU"

# step is only relevant for inserts derived from RE-digested input DNAs
if [ "$EXPECTING_ENDPOINT_RE_SITES" = "TRUE" ]; then

    # must always create a sample-level filtering site index
    echo "creating sample-level site lookup files"
    perl ${ACTION_DIR}/locate/create_index.pl | 
    $BGZIP $FILTERING_SITES_BGZ
    checkPipe

    echo "creating tabix index"
    $TABIX --sequence 1 --begin 2 --end 2 $FILTERING_SITES_BGZ
    checkPipe

    echo "done"

elif [ "$REJECTING_JUNCTION_RE_SITES" = "TRUE" ]; then

    # a genome-level filtering site index only needs to be created once per genome+enzyme combination
    if [ ! -f "$GENOME_FILTERING_SITES_BGZ" ]; then
        echo "creating genome-level site lookup files"
        mkdir -p $GENOME_ENZYME_DIR
        perl ${ACTION_DIR}/locate/create_index.pl | 
        $BGZIP $GENOME_FILTERING_SITES_BGZ
        checkPipe

        echo "creating tabix index"
        $TABIX --sequence 1 --begin 2 --end 2 $GENOME_FILTERING_SITES_BGZ
        checkPipe

        echo "done"
    else 
        echo "skipping genome-level site indexing, index already exists"
    fi

else 
    echo "skipping RE site indexing, not applicable to this library type"
fi

#!/bin/bash

# parse genome options into the requested genome

# function to display alert messages
show_conversion_message () {
    local GRC_GENOME=$1
    local UCSC_GENOME=$2
    echo
    echo "genome name ${GRC_GENOME} coerced to UCSC name ${UCSC_GENOME}"
    echo
}

# coerce GRC to UCSC genome names
# using UCSC names helps support the MDI genome browser
#   which dynamically reads genome metadata from the UCSC API
if [ "$GENOME" = "CHM13" ]; then
    show_conversion_message $GENOME hs1
    export GENOME=hs1
elif [ "$GENOME" = "GRCh38" ]; then
    show_conversion_message $GENOME hg38
    export GENOME=hg38
elif [ "$GENOME" = "GRCm39" ]; then
    show_conversion_message $GENOME mm39
    export GENOME=mm39
fi

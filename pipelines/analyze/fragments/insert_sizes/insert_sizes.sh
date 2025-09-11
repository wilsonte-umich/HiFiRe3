#!/bin/bash

# variables
INSERT_SIZES_DIR=${ACTION_DIR}/insert_sizes
if [[ "$CHECK_ENDPOINT_RE_MATCH" != "" && "$ENZYME_NAME" != "NA" ]]; then
    INSERT_SIZES_SCRIPT=${INSERT_SIZES_DIR}/insert_sizes.R
else 
    INSERT_SIZES_SCRIPT=${INSERT_SIZES_DIR}/non_RE/insert_sizes.R
fi

# analyze insert sizes
mkdir -p $PLOTS_DIR
Rscript ${INSERT_SIZES_SCRIPT}
checkPipe

echo "done"
echo

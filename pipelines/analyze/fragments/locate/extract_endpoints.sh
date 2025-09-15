#!/bin/bash

# extract observed sequence endpoints likely to match a RE cut site position
# group and count by unique endpoint

perl ${ACTION_DIR}/locate/extract_endpoints.pl | 
pigz -c --processes ${N_CPU} > ${OBSERVED_ENDPOINTS_FILE}
checkPipe

echo
echo "done"

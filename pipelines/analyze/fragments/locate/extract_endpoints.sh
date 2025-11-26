#!/bin/bash

# extract observed sequence endpoints likely to match a RE cut site position
# group and count by unique endpoint

# step is only relevant for reads expected to match RE site at their ends
if [ "$EXPECTING_ENDPOINT_RE_SITES" = "TRUE" ]; then

    perl ${ACTION_DIR}/locate/extract_endpoints.pl | 
    pigz -c --processes ${N_CPU} > ${OBSERVED_ENDPOINTS_FILE}
    checkPipe
    echo
    echo "done"
    
else 
    echo "skipping RE endpoint extraction, not applicable to this library type"
fi

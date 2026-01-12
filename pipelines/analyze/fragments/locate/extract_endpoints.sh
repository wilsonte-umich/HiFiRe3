#!/bin/bash

# extract observed sequence endpoints likely to match a RE cut site position
# group and count by unique endpoint

# step is only relevant for reads expected to match RE site at their ends
if [ "$EXPECTING_ENDPOINT_RE_SITES" = "TRUE" ] && [ "$CREATING_SAMPLE_SITE_FILES" = "TRUE" ]; then

    # perl ${ACTION_DIR}/locate/extract_endpoints.pl | 
    samtools view $NAME_BAM_FILE | 
    ${SUITE_BIN_DIR}/hf3_tools extract_endpoints | 
    pigz -c --processes ${N_CPU} > ${OBSERVED_ENDPOINTS_FILE}
    checkPipe
    echo "done"

elif [ "$EXPECTING_ENDPOINT_RE_SITES" = "TRUE" ]; then
    if [ "$SKIP_RFLP_DETECTION" = "1" ]; then
        echo "skipping RE endpoint extraction per user request, defaulting to in silico sites only"
    else 
        echo "skipping RE endpoint extraction, using prior site calls from --site-override-dir"
    fi
else 
    echo "skipping RE endpoint extraction, not applicable to this library type"
fi

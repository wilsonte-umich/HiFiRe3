#!/bin/bash

# build a table of observed endpoints intersected with in silico RE sites

# step is only relevant for reads expected to match RE site at their ends
if [ "$EXPECTING_ENDPOINT_RE_SITES" = "TRUE" ] && [ "$CREATING_SAMPLE_SITE_FILES" = "TRUE" ]; then

    Rscript ${ACTION_DIR}/locate/tabulate_endpoints.R
    checkPipe
    echo "done"

elif [ "$EXPECTING_ENDPOINT_RE_SITES" = "TRUE" ]; then
    if [ "$SKIP_RFLP_DETECTION" = "1" ]; then
        echo "skipping RE endpoint tabulation per user request, defaulting to in silico sites only"
    else 
        echo "skipping RE endpoint tabulation, using prior site calls from --site-override-dir"
    fi
else 
    echo "skipping RE endpoint tabulation, not applicable to this library type"
    # will index and use GENOME_SITES_GZ in place of FILTERING_SITES_FILE if REJECTING_JUNCTION_RE_SITES==TRUE
fi

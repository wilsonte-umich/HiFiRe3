#!/bin/bash

# build a table of observed endpoints intersected with in silico RE sites

# step is only relevant for reads expected to match RE site at their ends
if [ "$EXPECTING_ENDPOINT_RE_SITES" = "TRUE" ]; then

    Rscript ${ACTION_DIR}/locate/tabulate_endpoints.R
    checkPipe
    
else 
    echo "skipping RE endpoint tabulation, not applicable to this library type"
    # will index and use GENOME_SITES_GZ in place of FILTERING_SITES_FILE if REJECTING_JUNCTION_RE_SITES==TRUE
fi

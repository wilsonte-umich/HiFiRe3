#!/bin/bash

# build a table of observed endpoints intersected with in silico RE sites

# step is only relevant for reads expected to match RE site at their ends, when not overridden by options
if [ "$EXPECTING_ENDPOINT_RE_SITES" = "TRUE" ] && [ "$CREATING_SAMPLE_SITE_FILES" = "TRUE" ]; then
    Rscript ${ACTION_DIR}/locate/tabulate_endpoints.R | 
    awk 'BEGIN{printf "#"}{print}' | # add comment on header line
    bgzip -c > ${FILTERING_SITES_FILE}
    checkPipe
    tabix --sequence 1 --begin 2 --end 2 ${FILTERING_SITES_FILE}
    checkPipe
    echo "done"

# if expecting endpoint RE sites, user can override default behavior by skipping rflp detection...
elif [ "$EXPECTING_ENDPOINT_RE_SITES" = "TRUE" ]; then
    if [ "$SKIP_RFLP_DETECTION" = "1" ]; then
        echo "skipping RE endpoint tabulation per user request, defaulting to in silico sites only"
    
    # ... or by providing prior site calls from another pre-established file
    else 
        echo "skipping RE endpoint tabulation per user request, using sites listed in --site-override-file"
    fi

# otherwise index and use GENOME_SITES_BGZ==GENOME_FILTERING_SITES_FILE in place of FILTERING_SITES_FILE if REJECTING_JUNCTION_RE_SITES==TRUE
else 
    echo "skipping RE endpoint tabulation, not applicable to this library type"
fi

#!/bin/bash

# build a table of observed endpoints intersected with in silico RE sites
Rscript ${PIPELINE_SHARED}/tabulate_endpoints.R
checkPipe

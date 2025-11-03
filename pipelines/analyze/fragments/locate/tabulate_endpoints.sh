#!/bin/bash

# build a table of observed endpoints intersected with in silico RE sites
Rscript ${ACTION_DIR}/locate/tabulate_endpoints.R
checkPipe

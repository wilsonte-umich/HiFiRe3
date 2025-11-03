#!/bin/bash

# build a table of observed endpoints intersected with in silico RE sites
Rscript ${ACTION_DIR}/insert_sizes/analyze_insert_sizes.R
checkPipe

#!/bin/bash

# create BED files of all in silico RE sites for a set of REs
cat $GENOME_FASTA | 
perl ${MODULES_DIR}/genome/digest.pl
checkPipe

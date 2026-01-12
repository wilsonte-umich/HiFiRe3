#!/bin/bash

# reformat ONT unaligned bam files to compress them and to convert legacy files
runWorkflowStep 1 reformat $ACTION_DIR/reformat.sh

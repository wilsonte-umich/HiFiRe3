#!/bin/bash

# top-level script to coordinate read alignment to a haploid reference

# set working directory to READ_FILE_DIR to avoid too-long argument list with multiple read files
if [ "$READ_FILE_DIR" == "null" ]; then
    READ_FILE_DIR=${TASK_DIR}
fi
cd $READ_FILE_DIR

# set derivative environment variables and file paths
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/REs/set_re_vars.sh
source ${MODULES_DIR}/agFree/set_library_vars.sh
source $MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/align/set_read_file_vars.sh

# set sort and other parameters
source $MODULES_DIR/utilities/shell/create_temp_dir_small.sh
export SORT_RAM_PER_CPU_INT=$(((TOTAL_RAM_INT - 50000000000) / N_CPU))
export SAMTOOLS_CPU=$(( N_CPU - 1 )) # samtools uses 1 core for IO

# set the aligner
export ALIGNER=minimap2
export ALN_CPU=$(( N_CPU - 1 )) # minimap2 uses 1 core for IO
export ALIGN_SCRIPT=${ACTION_DIR}/align.sh # minimap2 implicitly handles admixed single and paired reads
export ALIGNER_OUTPUT_FILE=${NAME_BAM_FILE}
export MINIMAP2_INDEX=${MINIMAP2_INDEX_DIR}/${GENOME}.${ALIGNMENT_MODE}.mmi

# align read sequences as requested
runWorkflowStep 1 align ${PIPELINE_SHARED}/align.sh

# clean up
rm -r $TMP_DIR_WRK_SMALL

# initialize actions common to analyze actions

# set derivative environment variables and file paths
source $MODULES_DIR/genome/set_genome_vars.sh
# source $MODULES_DIR/align/set_alignment_vars.sh
# source $MODULES_DIR/REs/set_re_vars.sh
# source $MODULES_DIR/library/set_library_vars.sh
source $MODULES_DIR/analyze/set_analysis_vars.sh

# collect the working sample directories
if [ "$SAMPLE_DIRS" != "null" ]; then
    export SAMPLE_DIRS=$(echo $SAMPLE_DIRS | tr ',' ' ')
else 
    export SAMPLE_DIRS=$(ls -d $PROJECT_DIR/*/*.$VARIANT_FILE_SUFFIX | xargs -n1 dirname)
fi
export N_SAMPLES=$(echo $SAMPLE_DIRS | wc -w)
if [ $N_SAMPLES -lt 2 ]; then
    echo "ERROR: at least two sample directories are required for comparison"
    echo "check options --project-dir and --sample-dirs to be sure they point at `analyze $VARIANT_TYPE` output directories"
    exit 1
fi
export SAMPLE_FILES=$(echo "$SAMPLE_DIRS" | xargs -n 1 -I {} echo {}/*.$VARIANT_FILE_SUFFIX)

# set temporary directories
source $MODULES_DIR/utilities/shell/create_temp_dir_small.sh
source $MODULES_DIR/utilities/shell/create_temp_dir_shm.sh

# set other working parameters
export SORT_RAM_PER_CPU_INT=$(((TOTAL_RAM_INT - 50000000000) / N_CPU))
export SAMTOOLS_CPU=$((N_CPU - 1))
export ALN_CPU=$(( N_CPU - 1 )) # minimap2 uses 1 core for IO

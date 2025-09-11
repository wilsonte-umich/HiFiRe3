# initialize actions common to analyze actions

# set derivative environment variables and file paths
source $MODULES_DIR/genome/set_genome_vars.sh
source $MODULES_DIR/genotype/set_genotype_vars.sh
source $MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/REs/set_re_vars.sh
source $MODULES_DIR/agFree/set_agFree_vars.sh
source $MODULES_DIR/zoo/set_zoo_vars.sh
source $MODULES_DIR/analyze/set_analysis_vars.sh

# set temporary directories
source $MODULES_DIR/utilities/shell/create_temp_dir_small.sh
source $MODULES_DIR/utilities/shell/create_temp_dir_shm.sh

# set other working parameters
export SAMTOOLS_CPU=$((N_CPU - 1))

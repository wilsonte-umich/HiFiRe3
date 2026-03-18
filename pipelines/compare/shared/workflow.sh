# initialize actions common to compare actions

# set derivative environment variables and file paths
source $MODULES_DIR/rust/set_rust_vars.sh
source $MODULES_DIR/genome/set_genome_vars.sh
source $MODULES_DIR/align/set_alignment_vars.sh
source $MODULES_DIR/REs/set_re_vars.sh
source $MODULES_DIR/library/set_library_vars.sh
source $MODULES_DIR/analyze/set_analysis_vars.sh
source $MODULES_DIR/compare/set_compare_vars.sh

# set temporary directories
source $MODULES_DIR/utilities/shell/create_temp_dir_small.sh
source $MODULES_DIR/utilities/shell/create_temp_dir_shm.sh

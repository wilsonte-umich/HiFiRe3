# action:
#     set environment variables that identify and characterize input unaligned bam read files
# expects:
#     directory has been set upstream to $READ_FILE_DIR
#     sequence read file(s) in $READ_FILE_DIR in one of these patterns: (used in this precedence order)
#     if READ_FILE_PREFIX == NA
#         *.bam
#     else 
#         $READ_FILE_PREFIX*.bam
# usage:
#     source $MODULES_DIR/align/set_read_file_vars.sh

# set the sequence read input files
if [ "$READ_FILE_DIR" == "null" ]; then
    export READ_FILE_DIR=${TASK_DIR}/ubam
fi
if [ "$READ_FILE_PREFIX" == "null" ]; then
    export READ_FILE_PREFIX=""
fi

# find all matching read files within READ_FILE_DIR (default TASK_DIR)
export READ_FILES=`ls -1 $READ_FILE_DIR/$READ_FILE_PREFIX*.bam 2>/dev/null`
if [ "$READ_FILES" = "" ]; then
    echo "no valid input read files found in $READ_FILE_DIR"
    echo "check your values for options --read-file-dir and --read-file-prefix"
    exit 1
fi
export N_READ_FILES=`echo $READ_FILES | wc -w`

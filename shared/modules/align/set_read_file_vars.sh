# action:
#     apply default of TASK_DIR/ubam to READ_FILE_DIR if not set by user
#     set environment variables that identify and characterize input read files
# expects:
#     sequence read file(s) in $READ_FILE_DIR in one of these patterns: (used in this precedence order)
#     if READ_FILE_PREFIX == null
#         *.unaligned.bam
#         *.fastq.gz
#         *.fastq
#         *.bam
#     else 
#         $READ_FILE_PREFIX*.unaligned.bam
#         $READ_FILE_PREFIX*.fastq.gz
#         $READ_FILE_PREFIX*.fastq
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

# find all matching read files within READ_FILE_DIR
export READ_FILE_TYPE=unaligned.bam
export READ_1_FILES=`ls -1 $READ_FILE_DIR/$READ_FILE_PREFIX*.$READ_FILE_TYPE 2>/dev/null`
if [ "$READ_1_FILES" = "" ]; then
    export READ_FILE_TYPE=fastq.gz
    export READ_1_FILES=`ls -1 $READ_FILE_DIR/$READ_FILE_PREFIX*.$READ_FILE_TYPE 2>/dev/null`
fi
if [ "$READ_1_FILES" = "" ]; then
    export READ_FILE_TYPE=fastq
    export READ_1_FILES=`ls -1 $READ_FILE_DIR/$READ_FILE_PREFIX*.$READ_FILE_TYPE 2>/dev/null`
fi
if [ "$READ_1_FILES" = "" ]; then
    export READ_FILE_TYPE=bam
    export READ_1_FILES=`ls -1 $READ_FILE_DIR/$READ_FILE_PREFIX*.$READ_FILE_TYPE 2>/dev/null`
fi
if [ "$READ_1_FILES" = "" ]; then
    echo "no valid input read files found in $READ_FILE_DIR"
    echo "check your values for options --read-file-dir, --read-file-prefix, and --read-number-format"
    exit 1
fi
export N_READ_1_FILES=`echo $READ_1_FILES | wc -w`

# parse paired read files
export READ_2_FILES=""
export N_READ_2_FILES=0
if [ "$READ_PAIR_TYPE" = "paired" ]; then
    if [ $N_READ_1_FILES -gt 2 ]; then # with multiple files per read, a portion of the file name must identify the read number
        READ_MATCH=`echo $READ_NUMBER_FORMAT | sed 's/x/2/'`
        export READ_2_FILES=`ls $READ_1_FILES 2>/dev/null | grep $READ_MATCH | tr "\n" " "`
        READ_MATCH=`echo $READ_NUMBER_FORMAT | sed 's/x/1/'`
        export READ_1_FILES=`ls $READ_1_FILES 2>/dev/null | grep $READ_MATCH | tr "\n" " "`
    else # otherwise, if only two matching files, they can be named anything
        export READ_2_FILES=`echo $READ_1_FILES | cut -d " " -f2`
        export READ_1_FILES=`echo $READ_1_FILES | cut -d " " -f1`
    fi
    export N_READ_1_FILES=`echo $READ_1_FILES | wc -w`
    export N_READ_2_FILES=`echo $READ_2_FILES | wc -w`
    if [ "$N_READ_1_FILES" != "$N_READ_2_FILES" ]; then
        echo "different number of files found for paired reads 1 and 2"
        exit 1
    fi
fi

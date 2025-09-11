# action:
#     set environment variables that identify and characterize input read files
# expects:
#     directory has been set upstream to $READ_FILE_DIR
#     sequence read file(s) in $READ_FILE_DIR in one of these patterns: (used in this precedence order)
#     if READ_FILE_PREFIX == NA
#         *.fastq.gz
#         *.unaligned.bam
#         *.sra
#     else 
#         $READ_FILE_PREFIX*.fastq.gz
#         $READ_FILE_PREFIX*.unaligned.bam
#         $READ_FILE_PREFIX*.sra
#     if more than two files match, for paired reads, files must include read number identifiers as
#         *$READ_NUMBER_FORMAT*.fastq.gz 
#     e.g., for READ_NUMBER_FORMAT == _Rx_
#         *_R1_*.fastq.gz
#         *_R2_*.fastq.gz
# usage:
#     source $MODULES_DIR/source/set_read_file_vars.sh

# set the sequence read input files
if [ "$READ_FILE_PREFIX" == "null" ]; then
    export READ_FILE_PREFIX=""
fi

# find all matching read files within READ_FILE_DIR (default TASK_DIR)
export READ_FILE_TYPE=fastq.gz
export READ_1_FILES=`ls -1 $READ_FILE_PREFIX*.$READ_FILE_TYPE 2>/dev/null`
if [ "$READ_1_FILES" = "" ]; then
    export READ_FILE_TYPE=fastq
    export READ_1_FILES=`ls -1 $READ_FILE_PREFIX*.$READ_FILE_TYPE 2>/dev/null`
fi
if [ "$READ_1_FILES" = "" ]; then
    export READ_FILE_TYPE=unaligned.bam
    export READ_1_FILES=`ls -1 $READ_FILE_PREFIX*.$READ_FILE_TYPE 2>/dev/null`
fi
if [ "$READ_1_FILES" = "" ]; then # pattern typical of an ONT run with basecalling
    export READ_FILE_TYPE=unaligned.bam
    export READ_1_FILES=`ls -1 ubam/$READ_FILE_PREFIX*.$READ_FILE_TYPE 2>/dev/null`
fi
if [ "$READ_1_FILES" = "" ]; then
    export READ_FILE_TYPE=sra
    export READ_1_FILES=`ls -1 $READ_FILE_PREFIX*.$READ_FILE_TYPE 2>/dev/null`
fi
if [ "$READ_1_FILES" = "" ]; then
    echo "no valid input read files found in $PWD"
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

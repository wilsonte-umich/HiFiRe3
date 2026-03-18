#!/bin/bash

# run `dorado basecaller` on a cluster compute node, with optimized file transfer
# after dorado, run custom adapter trimming that is RE-site aware

# this newer version of the script uses POD5 auto-batching, to deal with progressively smaller pod5 file size as runs proceed
# input POD5 files will be grouped so that no two sequential batches will fill more than ~80% of POD5_BUFFER_SIZE
# which accounts for space utilization, i.e., that POD5_BUFFER_DIR must have free space matching:
#   two input batches of POD5 files
#   plus the smaller space required to store one batch's ubam files
# the process isn't perfect with now-larger POD5 files so ensure you have sufficient RAM if using RAM buffer

# parse option values
if [ "$DORADO_OPTIONS" == "null" ]; then DORADO_OPTIONS=""; fi
if [[ "$FORCE_BASECALLING" != "" && "$FORCE_BASECALLING" != "0" ]]; then FORCE_BASECALLING="true"; fi

# add modified base models
if [[ "$MODIFIED_BASE_MODELS" != "" && "$MODIFIED_BASE_MODELS" != "NA" && "$MODIFIED_BASE_MODELS" != "null" ]]; then 
    ONT_MODEL_COMPLEX="${ONT_MODEL_COMPLEX},${MODIFIED_BASE_MODELS}"
fi

# begin log report
echo "calling bases"
echo "  Dorado version:   "${DORADO_VERSION}
echo "  model(s):         "${ONT_MODEL_COMPLEX}
echo "  with options:     "${DORADO_OPTIONS}
echo "  pod5 buffer:      "${POD5_BUFFER_DIR}
echo "  buffer size:      "${POD5_BUFFER_SIZE}
echo "  input(s):         "${EXPANDED_INPUT_DIR}
echo "  output folder:    "${UBAM_DIR}

# prepare the cache and output directories
mkdir -p ${POD5_BUFFER_DIR}
mkdir -p ${UBAM_DIR}
rm -rf ${POD5_BUFFER_DIR}/*
BATCH_INDEX_FILE=$POD5_BUFFER_DIR/batches.txt

# functions for copy and calling a batch of pod5 files
do_batch_copy () {
    if [ "$COPY_OUT_I" != "" ]; then
        COPY_OUT_DIR=$POD5_BUFFER_DIR/${BATCH_PREFIX}_$COPY_OUT_I
        BATCH_OUTPUT_FILE2=$UBAM_DIR/${BATCH_PREFIX}_$COPY_OUT_I.unaligned.bam
        if [[ ! -f $BATCH_OUTPUT_FILE2 || "$FORCE_BASECALLING" == "true" ]]; then
            samtools view --with-header $COPY_OUT_DIR/out/*.bam | 
            # perl $MODULES_DIR/ONT/trim.pl | 
            ${HF3_TOOLS_BIN} trim_ont | 
            samtools view --bam --output $BATCH_OUTPUT_FILE2 -
            # mv -f $COPY_OUT_DIR/out/*.bam $UBAM_DIR
        fi
        rm -fr $COPY_OUT_DIR
    fi
    if [ "$COPY_IN_I" != "" ]; then
        BATCH_OUTPUT_FILE1=$UBAM_DIR/${BATCH_PREFIX}_$COPY_IN_I.unaligned.bam
        if [[ ! -f $BATCH_OUTPUT_FILE1 || "$FORCE_BASECALLING" == "true" ]]; then
            COPY_IN_DIR=$POD5_BUFFER_DIR/${BATCH_PREFIX}_$COPY_IN_I
            mkdir -p $COPY_IN_DIR/in
            mkdir -p $COPY_IN_DIR/out
            COPY_IN_BATCH_N=$(( COPY_IN_I + 1 ))
            echo "copying file batch $COPY_IN_BATCH_N"
            awk '$1=='$COPY_IN_BATCH_N $BATCH_INDEX_FILE
            BATCH_FILES=`awk '$1=='$COPY_IN_BATCH_N'{ print $3 }' $BATCH_INDEX_FILE`
            if [ "$IS_TAR_ARCHIVE" == "TRUE" ]; then
                tar -xf $TAR_ARCHIVE -C $COPY_IN_DIR/in --strip-components=1 $BATCH_FILES
            else 
                cp $BATCH_FILES $COPY_IN_DIR/in
            fi

        else 
            echo `basename $BATCH_OUTPUT_FILE1`" exists, set --force-basecalling to overwrite"
        fi
    fi
}
RUN_DORADO="$DORADO_EXECUTABLE basecaller --no-trim $DORADO_OPTIONS --models-directory ${DORADO_MODELS_DIR} $ONT_MODEL_COMPLEX" # we will control trimming downstream
run_batch_process () {
    PROCESS_DIR=$POD5_BUFFER_DIR/${BATCH_PREFIX}_$PROCESS_I
    BATCH_OUTPUT_FILE3=$UBAM_DIR/${BATCH_PREFIX}_$PROCESS_I.unaligned.bam
    TMP_OUPUT_FILE=$PROCESS_DIR/out/${BATCH_PREFIX}_$PROCESS_I.unaligned.bam
    if [[ ! -f $BATCH_OUTPUT_FILE3 || "$FORCE_BASECALLING" == "true" ]]; then
        $RUN_DORADO $PROCESS_DIR/in > $TMP_OUPUT_FILE
    fi
    rm -fr $PROCESS_DIR/in
}

# process one pod5 input directory at a time; output is merged into a single output directory
WORKING_I=0
for WORKING_INPUT_DIR in ${EXPANDED_INPUT_DIR}; do

echo
echo "processing $WORKING_INPUT_DIR"
WORKING_I=$((WORKING_I + 1))
BATCH_PREFIX="batch_$WORKING_I"

# initialize pod5 sources
if [[ $WORKING_INPUT_DIR == *.tar ]]; then
    IS_TAR_ARCHIVE=TRUE
fi
if [ "$IS_TAR_ARCHIVE" == "TRUE" ]; then
    TAR_ARCHIVE=`basename $WORKING_INPUT_DIR`
    TAR_INDEX="$TAR_ARCHIVE.contents"
    WORKING_INPUT_DIR=`dirname $WORKING_INPUT_DIR`
    cd ${WORKING_INPUT_DIR}
    if [ ! -f $TAR_INDEX ]; then
        echo "listing tar archive contents" 
        tar -tvf $TAR_ARCHIVE | grep -P '\.pod5$' > $TAR_INDEX # -v is like ls -l
    fi
    N_POD5_BATCHES=`cat $TAR_INDEX | perl $ACTION_DIR/autobatch_pod5.pl 2 $BATCH_INDEX_FILE`
else 
    cd ${WORKING_INPUT_DIR}
    N_POD5_BATCHES=`ls -l *.pod5 | perl $ACTION_DIR/autobatch_pod5.pl 4 $BATCH_INDEX_FILE`
fi
echo "processing files in $N_POD5_BATCHES batches"

# run basecaller one pod5 batch at a time, working from /dev/shm
# once the loop is running, one batch is copying while another batch is basecalling
COPY_IN_I=0
COPY_OUT_I=""
echo "waiting for batch copy $COPY_IN_I"
do_batch_copy
PROCESS_I=$COPY_IN_I
for ((COPY_IN_I=1; COPY_IN_I < $N_POD5_BATCHES; COPY_IN_I+=1)); do 
    do_batch_copy &
    COPY_PID=$!
    run_batch_process
    echo "waiting for batch copy $COPY_IN_I"
    wait $COPY_PID
    COPY_OUT_I=$PROCESS_I
    PROCESS_I=$COPY_IN_I
done
COPY_IN_I=""
do_batch_copy &
COPY_PID=$!
run_batch_process
echo "finishing final file copy"
wait $COPY_PID
COPY_OUT_I=$PROCESS_I
do_batch_copy

# end input directory loop
done


# VERSION 0.5.3
# important changes: 

# added --trim/--no-trim options to basecaller only (not duplex)
# however, duplex reads will have the adapters trimmed off anyway per:
#   https://github.com/nanoporetech/dorado/issues/509
#   https://github.com/nanoporetech/dorado/issues/510
# also note that: chimeric read splitting is enabled for both duplex and simplex basecalling by default
# however, some chimeric reads with unsplittable junctions escape duplex detection

# automated model selection, although in general we don't intend to use this feature:
#    model {fast,hac,sup}@v{version} for automatic model selection including modbases, or path to existing model directory
#    (fast|hac|sup)[@(version|latest)][,modification[@(version|latest)]][,...]


# VERSION 0.7.1
# note addition of --mm2-preset, to use the new minimap2 preset for >Q20 long reads, lr:hq

# $ ./dorado basecaller --help
# [2024-06-17 16:38:07.836] [info] Running: "basecaller" "--help"
# Usage: dorado [-h] [--device VAR] [--read-ids VAR] [--resume-from VAR] [--max-reads VAR] [--min-qscore VAR] [--batchsize VAR] [--chunksize VAR] [--overlap VAR] [--recursive] [--modified-bases VAR...] [--modified-bases-models VAR] [--modified-bases-threshold VAR] [--emit-fastq] [--emit-sam] [--emit-moves] [--reference VAR] [--kit-name VAR] [--barcode-both-ends] [--no-trim] [--trim VAR] [--sample-sheet VAR] [--barcode-arrangement VAR] [--barcode-sequences VAR] [--primer-sequences VAR] [--estimate-poly-a] [--poly-a-config VAR] [-k VAR] [-w VAR] [-I VAR] [--secondary VAR] [-N VAR] [-Y] [--bandwidth VAR] [--junc-bed VAR] [--mm2-preset VAR] model data

# Positional arguments:
#   model                         model selection {fast,hac,sup}@v{version} for automatic model selection including modbases, or path to existing model directory 
#   data                          the data directory or file (POD5/FAST5 format). 

# Optional arguments:
#   -h, --help                    shows help message and exits 
#   -v, --verbose             
#   -x, --device                  device string in format "cuda:0,...,N", "cuda:all", "metal", "cpu" etc.. [default: "cuda:all"]
#   -l, --read-ids                A file with a newline-delimited list of reads to basecall. If not provided, all reads will be basecalled [default: ""]
#   --resume-from                 Resume basecalling from the given HTS file. Fully written read records are not processed again. [default: ""]
#   -n, --max-reads               [default: 0]
#   --min-qscore                  Discard reads with mean Q-score below this threshold. [default: 0]
#   -b, --batchsize               if 0 an optimal batchsize will be selected. batchsizes are rounded to the closest multiple of 64. [default: 0]
#   -c, --chunksize               [default: 10000]
#   -o, --overlap                 [default: 500]
#   -r, --recursive               Recursively scan through directories to load FAST5 and POD5 files 
#   --modified-bases              [nargs: 1 or more] 
#   --modified-bases-models       a comma separated list of modified base models [default: ""]
#   --modified-bases-threshold    the minimum predicted methylation probability for a modified base to be emitted in an all-context model, [0, 1] [default: 0.05]
#   --emit-fastq                  Output in fastq format. 
#   --emit-sam                    Output in SAM format. 
#   --emit-moves              
#   --reference                   Path to reference for alignment. [default: ""]
#   --kit-name                    Enable barcoding with the provided kit name. Choose from: EXP-NBD103 EXP-NBD104 EXP-NBD114 EXP-NBD196 EXP-PBC001 EXP-PBC096 SQK-16S024 SQK-16S114-24 SQK-LWB001 SQK-MLK111-96-XL SQK-MLK114-96-XL SQK-NBD111-24 SQK-NBD111-96 SQK-NBD114-24 SQK-NBD114-96 SQK-PBK004 SQK-PCB109 SQK-PCB110 SQK-PCB111-24 SQK-PCB114-24 SQK-RAB201 SQK-RAB204 SQK-RBK001 SQK-RBK004 SQK-RBK110-96 SQK-RBK111-24 SQK-RBK111-96 SQK-RBK114-24 SQK-RBK114-96 SQK-RLB001 SQK-RPB004 SQK-RPB114-24 TWIST-16-UDI TWIST-96A-UDI VSK-PTC001 VSK-VMK001 VSK-VMK004 VSK-VPS001. [default: ""]
#   --barcode-both-ends           Require both ends of a read to be barcoded for a double ended barcode. 
#   --no-trim                     Skip trimming of barcodes, adapters, and primers. If option is not chosen, trimming of all three is enabled. 
#   --trim                        Specify what to trim. Options are 'none', 'all', 'adapters', and 'primers'. Default behaviour is to trim all detected adapters, primers, or barcodes. Choose 'adapters' to just trim adapters. The 'primers' choice will trim adapters and primers, but not barcodes. The 'none' choice is equivelent to using --no-trim. Note that this only applies to DNA. RNA adapters are always trimmed. [default: ""]
#   --sample-sheet                Path to the sample sheet to use. [default: ""]
#   --barcode-arrangement         Path to file with custom barcode arrangement. 
#   --barcode-sequences           Path to file with custom barcode sequences. 
#   --primer-sequences            Path to file with custom primer sequences. [default: <not representable>]
#   --estimate-poly-a             Estimate poly-A/T tail lengths (beta feature). Primarily meant for cDNA and dRNA use cases. 
#   --poly-a-config               Configuration file for PolyA estimation to change default behaviours [default: ""]
#   -k                            minimap2 k-mer size for alignment (maximum 28). 
#   -w                            minimap2 minimizer window size for alignment. 
#   -I                            minimap2 index batch size. 
#   --secondary                   minimap2 outputs secondary alignments 
#   -N                            minimap2 retains at most INT secondary alignments 
#   -Y                            minimap2 uses soft clipping for supplementary alignments 
#   --bandwidth                   minimap2 chaining/alignment bandwidth and optionally long-join bandwidth specified as NUM,[NUM] 
#   --junc-bed                    Optional file with gene annotations in the BED12 format (aka 12-column BED), or intron positions in 5-column BED. With this option, minimap2 prefers splicing in annotations. 
#   --mm2-preset                  minimap2 preset for indexing and mapping. Alias for the -x option in minimap2. [default: "lr:hq"]

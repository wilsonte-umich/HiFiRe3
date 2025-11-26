# actions:
#     align read sequences
# expects:
#     source $MODULES_DIR/genome/set_genome_vars.sh
#     source $MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/align/set_read_file_vars.sh
#     source $MODULES_DIR/library/set_library_vars.sh
# optional for all read types:
#     $FORCE_ALIGNMENT  [default: don't overwrite NAME_BAM_FILE]
# input:
#     $READ_1_FILES (and $READ_2_FILES if paired)
#     if ONT, reads subjected to trimming by modules/ONT/trim.pl (NOT by Dorado)
# output:
#     NAME_BAM_FILE
# where bam:
#     for ONT, includes basecalling tags (missing for other platforms)
#       ch:i: = ONT channel
#       tl:Z: = ONT adapter trim lengths in format `<5' trim>,<3' trim>`
#     includes cs:Z: tag from minimap2
#     includes custom alignment tags:
#       fm:Z: = read merge status (for paired reads)
#       hv:i: = bit-encoded flag of alignment and read variant status

#------------------------------------------------------------------
# set the product alignment file; abort silently if exists and not forced
#------------------------------------------------------------------
if [[ "$FORCE_ALIGNMENT" != "" && "$FORCE_ALIGNMENT" != "0" && "$FORCE_ALIGNMENT" != "false" && -e $NAME_BAM_FILE ]]; then
    echo "forcing overwrite of alignment file: "`basename $NAME_BAM_FILE`
    rm -f $NAME_BAM_FILE
fi
if [ -e $NAME_BAM_FILE ]; then
    echo "alignment file already exists: "`basename $NAME_BAM_FILE`
    echo "manually remove it or set --force-alignment to allow alignment to run again"
else

#------------------------------------------------------------------
# provide log feedback
#------------------------------------------------------------------
echo "aligning read sequences to $GENOME using minimap2"
echo "  platform: $SEQUENCING_PLATFORM"
echo "  library type: $LIBRARY_TYPE"
echo "  aln mode: $ALIGNMENT_MODE" 
echo "  bandwidth: $BANDWIDTH" 
echo "  reference genome: $GENOME_FASTA"
echo "  read files:"
echo "$READ_1_FILES" | perl -ne 'print "    ".join("\n    ", split(" ", $_)),"\n"'
if [ "$READ_PAIR_TYPE" = "paired" ]; then
    echo "$READ_2_FILES" | perl -ne 'print "    ".join("\n    ", split(" ", $_)),"\n"'
fi
echo "  output:"
echo "    $NAME_BAM_FILE"
echo

#------------------------------------------------------------------
# as needed, create and save the appropriate minimap2 genome alignment index
#------------------------------------------------------------------
export GENOME_FASTA_WRK=${GENOME_FASTA}
export ALIGNMENT_MODE_WRK=${ALIGNMENT_MODE}
export CREATE_MM2_INDEX_MESSAGE="continuing with alignment"
source ${MODULES_DIR}/align/create_mm2_index.sh # sets variable ${MINIMAP2_INDEX_WRK}

#------------------------------------------------------------------
# platform-specific fastp pre-alignment read processing with:
#       adapter trimming
#       read merging
#       read-level quality filtering
#   large numbers of threads do not improve fastp speed (its already fast)
#   dup evaluation is memory intensive and not needed
#   as desired, fastp prefers read1 bases over read2 bases when merging
#------------------------------------------------------------------
export FASTP_STEP_COMMAND=cat
export MERGE_STEP_COMMAND=cat
if [ "${RUN_PREALIGNMENT_FASTP}" != "" ]; then
    #   use buffered fastp until its memory leak is fixed: https://github.com/OpenGene/fastp/issues/392
    #   the leak was supposedly fixed in fastp v0.24 but memory climbing still occurs
    #   continue using buffered_fastp.pl as it adds little overhead since aligner is the slow step
    export FASTP_STEP_COMMAND="perl $ACTION_DIR/buffered_fastp.pl"
    export MERGE_STEP_COMMAND="perl $ACTION_DIR/adjust_merge_tags.pl"
    export FASTP_BUFFER=${TMP_DIR_WRK_SMALL}/FASTP_BUFFER.fastq
    export FASTP_BUFFER_N_READS=1000000
    export FASTP_COMMAND="fastp \
        --thread 3 --stdin --stdout --dont_eval_duplication \
        --length_required $PLATFORM_MIN_INSERT_SIZE \
        $FASTP_PAIRED_END_OPTIONS $FASTP_POLY_X_OPTIONS $FASTP_ADAPTER_OPTIONS \
        --qualified_quality_phred $QUALIFIED_QUALITY_PHRED --unqualified_percent_limit $UNQUALIFIED_PERCENT_LIMIT \
        --average_qual $AVERAGE_QUAL \
        --html $FASTP_LOG_PREFIX.html --json $FASTP_LOG_PREFIX.json \
        --report_title \"$DATA_NAME\" 2>$FASTP_LOG_PREFIX.txt"
    # --disable_quality_filtering
fi

#------------------------------------------------------------------
# process reads and align to genome
#------------------------------------------------------------------

# parse reads from ubam to fastq, tallying various metadata along the way
perl ${ACTION_DIR}/align/prepare_fastq.pl |

# if needed for short reads, use fastp for one-pass adapter trimming, quality filtering, and read merging
$FASTP_STEP_COMMAND | 

# tweak the way read pair merge status is reported
# report summary of fastp actions
$MERGE_STEP_COMMAND | 

# align to genome
#   soft-clip supplementary
#   include cs:Z: tag
#   copy FASTQ name extensions to SAM tags
minimap2 \
    -ax ${ALIGNMENT_MODE} \
    -t ${ALN_CPU} \
    -y -Y --secondary=no --cs \
    -r ${BANDWIDTH} \
    ${MINIMAP2_INDEX_WRK} - 2>${MINIMAP_LOG_FILE} |

# process reads to determine if variants are present according to the reference genome
perl ${ACTION_DIR}/align/flag_variant_events.pl |
samtools view -b -@ ${SAMTOOLS_CPU} > ${NAME_BAM_FILE}
checkPipe

echo "done"

fi

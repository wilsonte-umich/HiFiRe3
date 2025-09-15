# actions:
#     align read sequences
# expects:
#     source $MODULES_DIR/genome/set_genome_vars.sh
#     source $MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/source/set_read_file_vars.sh
#     source $MODULES_DIR/library/set_library_vars.sh
# optional for all read types:
#     $FORCE_ALIGNMENT  [default: don't overwrite ALIGNER_OUTPUT_FILE]
# input:
#     $READ_FILES
#     if ONT, reads subjected to trimming by modules/ONT/trim.pl (NOT by Dorado)
# output:
#     NAME_BAM_FILE
# where bam:
#     for ONT, includes basecalling tags (missing for PacBio)
#       ch:i: = ONT channel
#       tl:Z: = ONT adapter trim lengths in format `<5' trim>,<3' trim>`
#     includes cs:Z: tag from minimap2
#     includes custom alignment tags:
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
echo "  bandwidth: $BANDWIDTH" 
echo "  aln mode: $ALIGNMENT_MODE" 
echo "  reference genome: $GENOME_FASTA"
echo "  read files:"
echo "$READ_FILES" | perl -ne 'print "    ".join("\n    ", split(" ", $_)),"\n"'
echo "  output: $NAME_BAM_FILE"
echo

#------------------------------------------------------------------
# as needed, create and save the appropriate minimap2 genome alignment index
#------------------------------------------------------------------
export GENOME_FASTA_WRK=${GENOME_FASTA}
export ALIGNMENT_MODE_WRK=${ALIGNMENT_MODE}
source ${MODULES_DIR}/align/create_mm2_index.sh # sets variable ${MINIMAP2_INDEX_WRK}

#------------------------------------------------------------------
# process reads and align to genome
#------------------------------------------------------------------

# parse reads from ubam to fastq, tallying various metadata along the way
perl ${ACTION_DIR}/align/prepare_fastq.pl |

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

# action:
#   execute haploid reference alignment in a stream using minimap2
# input:
#   (interleaved) FASTQ stream on STDIN
# output:
#   NAME_BAM_FILE

# align to genome
#   soft-clip supplementary
#   include cs:Z: tag
minimap2 \
    -ax $ALIGNMENT_MODE \
    -t $ALN_CPU \
    -Y --secondary=no --cs \
    -r $BANDWIDTH \
    $MINIMAP2_INDEX - 2>$MINIMAP_LOG_FILE |

# process reads to determine if variants are present according to the reference genome
perl $ACTION_DIR/../shared/flag_variant_events.pl |
samtools view -b --threads $ALN_CPU --keep-tag NM,MD,AS,SA,cs,de,xf,xh | # --keep-tag does not change tag order
slurp -s 100M -o $NAME_BAM_FILE

# prepare for haploid reference alignment

# ------------------------------------------------------------------
# as needed, create and save the appropriate minimap2 alignment index
# ------------------------------------------------------------------
MINIMAP2_INDEX_LOCK=$MINIMAP2_INDEX.lock
while [ -f $MINIMAP2_INDEX_LOCK ]; do # prevent two jobs from creating the same index at the same time
    echo "waiting for minimap2 index to be created"
    sleep 60
done
if [ ! -f "$MINIMAP2_INDEX" ]; then
    mkdir -p $MINIMAP2_INDEX_DIR
    touch $MINIMAP2_INDEX_LOCK
    echo
    echo "creating minimap2 index for:"
    echo "    genome: "$GENOME
    echo "    fasta:  "`basename $GENOME_FASTA`
    echo "    mode:   "$ALIGNMENT_MODE
    minimap2 -x $ALIGNMENT_MODE -t 3 -d $MINIMAP2_INDEX $GENOME_FASTA
    checkPipe
    rm $MINIMAP2_INDEX_LOCK
fi

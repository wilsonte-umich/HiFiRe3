# actions:
#   combine SNV-prepared 5'-most alignments into read pileup by genome coordinate
# input:
#   $SNV_ALNS_PREFIX.<zeroPaddedChromIndex1>.txt.gz
# output:
#   tabix-indexed genome pileup of va-encoded alignments
#   tabix-indexed list of (sub)clonal SNVs and indels, optionally with expected genotypes

# process pileups
rm -f ${SNV_GENOME_PILEUP_PREFIX}.*
rm -f ${SNV_SUMMARY_TABLE_PREFIX}.*
perl ${ACTION_DIR}/pileup/pileup.pl
checkPipe

# combine and index the genome-level pileups and other files over all chrom threads
echo "concatenating genome-level pileups and SNV tables"
concatenate_chrom_files () {
    PREFIX=$1
    OUT_FILE=$2
    TYPE=$3
    echo "  "${TYPE}
    zcat ${PREFIX}.*.gz | 
    bgzip --threads ${N_CPU} --force --output ${OUT_FILE}
    checkPipe
    tabix --threads ${N_CPU} -p bed ${OUT_FILE}
    checkPipe
}
concatenate_chrom_files ${SNV_GENOME_PILEUP_PREFIX} ${SNV_GENOME_PILEUP} genome_pileup
concatenate_chrom_files ${SNV_SUMMARY_TABLE_PREFIX} ${SNV_SUMMARY_TABLE} snv_summary

echo "done"

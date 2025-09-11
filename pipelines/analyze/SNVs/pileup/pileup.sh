# actions:
#   combine SNV-prepared 5'-most alignments into two kinds of read pileups:
#       by genome coordinate
#       by position along:
#           read 1
#           read 2 (when applicable for paired end)
# input:
#   $SNV_ALNS_PREFIX.<zeroPaddedChromIndex1>.<strandIndex0>.txt.gz
# output:
#   tabix-indexed genome pileup of va-encoded alignments
#   tabix-indexed list of (sub)clonal SNVs and indels, optionally with expected genotypes
#   flat file of read pileups by read position relative to 5'-most end

# process pileups
rm -f ${SNV_GENOME_PILEUP_PREFIX}.*
rm -f ${SNV_SUMMARY_TABLE_PREFIX}.*
rm -f ${SNV_READ_PILEUP_PREFIX}.*
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

# aggregate the read-level pileups over all chrom threads
echo "concatenating read-level pileups"
concatenate_read_pileup(){
    PREFIX=$1
    OUT_FILE=$2
    TYPE=$3
    echo "  "${TYPE}
    (
        echo -e "readN\treadPos1\trefBases\taltBases\tzygosity\tgenotypeZygosity\tgenotypeCoverage\tnObserved"; 
        zcat ${PREFIX}.*.gz
    ) | gzip -c > ${OUT_FILE}
    checkPipe
}
concatenate_read_pileup ${SNV_READ_PILEUP_PREFIX} ${SNV_READ_PILEUP} read_pileup

echo "done"

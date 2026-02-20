# actions:
#   parse dd and cs tags to call SNV/indels and create a read pileup
#   see Rust crate `hf3_tools` for details
# input:
#   $NAME_BAM_FILES from one or more samples with hf3 tags
# output:
#   $SNV_PILEUP_BGZ
#   $SNV_VARIANTS_BGZ

# working variables
export INDEX_FILE_PREFIX_WRK=${TMP_FILE_PREFIX_SMALL}.index_snv
rm -f $INDEX_FILE_PREFIX_WRK.*.bed.bgz
rm -f $INDEX_FILE_PREFIX_WRK.*.txt.bgz
BGZIP="bgzip --threads $N_CPU --stdout"
TABIX="tabix --threads $N_CPU"

# split name-sorted BAM alignments by chromosome
# only on-target reads are retained for SNV calling
# this represents a first sort action and supports downstream parallelization by chrom
${SUITE_BIN_DIR}/hf3_tools split_bam_by_chrom_snv
checkPipe

# call SNVs and indels per chromosome
${SUITE_BIN_DIR}/hf3_tools analyze_snvs
checkPipe

# concatenate and index the pileups
echo "concatenating and indexing pileup files"
zcat $INDEX_FILE_PREFIX_WRK.chr*.all_reads.pileup.bed.bgz | 
$BGZIP > $SNV_ALL_READS_PILEUP_BGZ
checkPipe
$TABIX -p bed $SNV_ALL_READS_PILEUP_BGZ
checkPipe

zcat $INDEX_FILE_PREFIX_WRK.chr*.error_corrected.pileup.bed.bgz | 
$BGZIP > $SNV_ERROR_CORRECTED_PILEUP_BGZ
checkPipe
$TABIX -p bed $SNV_ERROR_CORRECTED_PILEUP_BGZ
checkPipe

# concatenate and index the allowed variant lists
echo "concatenating and indexing allowed variants files"
zcat $INDEX_FILE_PREFIX_WRK.chr*.all_reads.snv_indel.txt.bgz | 
$BGZIP > $SNV_ALL_READS_VARIANTS_BGZ
checkPipe
$TABIX --sequence 1 --begin 2 --end 2 $SNV_ALL_READS_VARIANTS_BGZ
checkPipe

zcat $INDEX_FILE_PREFIX_WRK.chr*.error_corrected.snv_indel.txt.bgz | 
$BGZIP > $SNV_ERROR_CORRECTED_VARIANTS_BGZ
checkPipe
$TABIX --sequence 1 --begin 2 --end 2 $SNV_ERROR_CORRECTED_VARIANTS_BGZ
checkPipe

echo
echo "done"

# actions:
#   call SNV/indels in PacBio by-strand sequencing
#   see Rust crate `hf3_tools` for details
# input:
#   $NAME_BAM_FILES from one or more samples with hf3 tags

# working variables
export INDEX_FILE_PREFIX_WRK=${TMP_FILE_PREFIX_SMALL}.index_snv
rm -f $INDEX_FILE_PREFIX_WRK.*.bgz
BGZIP="bgzip --threads $N_CPU --stdout"
TABIX="tabix --threads $N_CPU"

# split name-sorted BAM alignments by chromosome
# only on-target reads are retained for SNV calling
# this represents a first sort action and supports downstream parallelization by chrom
${HF3_TOOLS_BIN} split_bam_by_chrom_snv
checkPipe

# call SNVs and indels per chromosome
${HF3_TOOLS_BIN} analyze_snvs
checkPipe
echo

# concatenate and index the pileups
echo "concatenating and indexing pileup files"
# zcat $INDEX_FILE_PREFIX_WRK.chr*.all_reads.pileup.bed.bgz | 
# $BGZIP > $SNV_ALL_READS_PILEUP_BGZ
# checkPipe
# $TABIX -p bed $SNV_ALL_READS_PILEUP_BGZ
# checkPipe

zcat $INDEX_FILE_PREFIX_WRK.chr*.error_corrected.pileup.bed.bgz | 
$BGZIP > $SNV_ERROR_CORRECTED_PILEUP_BGZ
checkPipe
$TABIX -p bed $SNV_ERROR_CORRECTED_PILEUP_BGZ
checkPipe

# concatenate and index the allowed variant lists
echo "concatenating and indexing allowed variants files"
# zcat $INDEX_FILE_PREFIX_WRK.chr*.all_reads.snv_indel.txt.bgz | 
# $BGZIP > $SNV_ALL_READS_VARIANTS_BGZ
# checkPipe
# $TABIX --sequence 1 --begin 2 --end 2 $SNV_ALL_READS_VARIANTS_BGZ
# checkPipe

zcat $INDEX_FILE_PREFIX_WRK.chr*.error_corrected.snv_indel.txt.bgz | 
$BGZIP > $SNV_ERROR_CORRECTED_VARIANTS_BGZ
checkPipe
$TABIX --sequence 1 --begin 2 --end 2 $SNV_ERROR_CORRECTED_VARIANTS_BGZ
checkPipe

# concatenate and index the read encodings
echo "concatenating and indexing read encodings"
# zcat $INDEX_FILE_PREFIX_WRK.chr*.all_reads.read_encodings.bed.bgz | 
# $BGZIP > $SNV_ALL_READS_ENCODINGS_BGZ
# checkPipe
# $TABIX -p bed $SNV_ALL_READS_ENCODINGS_BGZ
# checkPipe

zcat $INDEX_FILE_PREFIX_WRK.chr*.error_corrected.read_encodings.bed.bgz | 
$BGZIP > $SNV_ERROR_CORRECTED_ENCODINGS_BGZ
checkPipe
$TABIX -p bed $SNV_ERROR_CORRECTED_ENCODINGS_BGZ
checkPipe

echo
echo "done"

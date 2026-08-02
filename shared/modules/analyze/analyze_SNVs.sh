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

# concatenate and index the allowed variant lists
echo "concatenating and indexing allowed variants file"
zcat $INDEX_FILE_PREFIX_WRK.chr*.snv_indel.variants.txt.bgz | 
$BGZIP > $SNV_VARIANTS_BGZ
checkPipe
$TABIX --sequence 1 --begin 2 --end 2 --zero-based $SNV_VARIANTS_BGZ
checkPipe

# concatenate and index the variant read list
echo "concatenating and indexing variant reads file"
zcat $INDEX_FILE_PREFIX_WRK.chr*.snv_indel.variant_reads.txt.bgz |
sort --parallel $N_CPU -S 4G -k1,1 -k2,2n -k3,3n | 
$BGZIP > $SNV_VARIANT_READS_BGZ
checkPipe
$TABIX --sequence 1 --begin 3 --end 4 --zero-based $SNV_VARIANT_READS_BGZ
checkPipe

# concatenate and index the read encodings
echo "concatenating and indexing fragment reads_on_ref encodings"
zcat $INDEX_FILE_PREFIX_WRK.chr*.fragments.on_reference.bed.bgz | 
$BGZIP > $SNV_FRAGMENTS_ON_REFERENCE_BGZ
checkPipe
$TABIX -p bed $SNV_FRAGMENTS_ON_REFERENCE_BGZ
checkPipe

echo "concatenating and indexing fragment reads_on_hap encodings"
zcat $INDEX_FILE_PREFIX_WRK.chr*.fragments.on_haplotype.bed.bgz | 
$BGZIP > $SNV_FRAGMENTS_ON_HAPLOTYPE_BGZ
checkPipe
$TABIX -p bed $SNV_FRAGMENTS_ON_HAPLOTYPE_BGZ
checkPipe

echo
echo "done"

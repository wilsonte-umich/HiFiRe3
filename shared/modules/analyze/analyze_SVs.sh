# actions:
#   extract all junctions and alignments from reads
#   fuzzy-match junctions to each other to yield final junction calls
#   see Rust crate `hf3_tools` for details
# input:
#   $NAME_BAM_FILES from one or more samples with hf3 tags
# output:
#   $SV_ALIGNMENTS_FILE
#   $SV_READ_PATHS_FILE
#   $SV_FINAL_JUNCTIONS_FILE 1 and 2

# working variables
export INDEX_FILE_PREFIX_WRK=${TMP_FILE_PREFIX_SMALL}.index_sv
rm -f $INDEX_FILE_PREFIX_WRK.*.bam
rm -f $INDEX_FILE_PREFIX_WRK.*.txt.gz
rm -f $INDEX_FILE_PREFIX_WRK.*.txt.bgz
TABIX="tabix --threads $N_CPU"

# split name-sorted BAM by the chromosome of the first alignment
# only on-target reads are retained for SV calling
# this represents a first sort action and supports downstream parallelization by chrom
${SUITE_BIN_DIR}/hf3_tools split_bam_by_chrom_sv
checkPipe

# index and parse unique fragment alignments and junctions
# fuzzy-match junctions to yield final SV calls
${SUITE_BIN_DIR}/hf3_tools analyze_svs
checkPipe
rm -f $INDEX_FILE_PREFIX_WRK.first_alns.chr*.gz

# index output files for app
echo "indexing read paths" # indexed by QNAME and read_len as a convenient way of indexed read retrieval
$TABIX --sequence 1 --begin 2 --end 2 --zero-based $SV_READ_PATHS_FILE
checkPipe

echo "indexing alignment map" # alignment map and coverage map do not have a column header
$TABIX --sequence 1 --begin 2 --end 2 $SV_ALIGNMENTS_FILE # NOT --begin 2 --end 3 because column 3 can be less than column 2 
checkPipe
echo "indexing coverage map"
$TABIX -p bed $SV_COVERAGE_FILE
checkPipe

echo "indexing final junctions by breakpoint 1"
$TABIX --sequence 1 --begin 2 --end 2 $SV_FINAL_JUNCTIONS_FILE_1
checkPipe
echo "indexing final junctions by breakpoint 2"
$TABIX --sequence 4 --begin 5 --end 5 $SV_FINAL_JUNCTIONS_FILE_2
checkPipe

echo
echo "done"

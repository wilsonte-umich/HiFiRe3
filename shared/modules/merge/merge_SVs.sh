# actions:
#   extract all previously called final junctions
#   fuzzy-match junctions to each other to yield merged final junction calls
#   see Rust crate `hf3_tools` for details
# input:
#   $FINAL_JUNCTION_FILES from one or more samples 
# output:
#   $SV_FINAL_JUNCTIONS_FILE 1 and 2
#   etc.

# working variables
TABIX="tabix --threads $N_CPU"

# fuzzy-match final junctions to yield final SV calls
${HF3_TOOLS_BIN} merge_svs
checkPipe

# index output files for app
echo "indexing read paths" # indexed by QNAME and read_len as a convenient way of indexed read retrieval
$TABIX --sequence 1 --begin 2 --end 2 --zero-based $SV_READ_PATHS_FILE
checkPipe

echo "indexing alignment map" # alignment map and coverage map do not have a column header
$TABIX --sequence 1 --begin 2 --end 2 $SV_ALIGNMENTS_FILE # NOT --begin 2 --end 3 because column 3 can be less than column 2 
checkPipe
# echo "indexing coverage map"
# $TABIX -p bed $SV_COVERAGE_FILE
# checkPipe

echo "indexing final junctions by breakpoint 1"
$TABIX --sequence 1 --begin 2 --end 2 $SV_FINAL_JUNCTIONS_FILE_1
checkPipe
echo "indexing final junctions by breakpoint 2"
$TABIX --sequence 4 --begin 5 --end 5 $SV_FINAL_JUNCTIONS_FILE_2
checkPipe

echo
echo "done"

# actions:
#   remove SEQ, QUAL and CIGAR from final junctions and read paths files
#   mask with *
#   the resulting simplified masked files are used for data packages for faster app loading
# input:
#   $SV_FINAL_JUNCTIONS_FILE_1
#   $SV_FINAL_JUNCTIONS_FILE_2

# output:
#   $SV_MASKED_JUNCTIONS_FILE_1
#   $SV_MASKED_JUNCTIONS_FILE_2

# working variables
BGZIP="bgzip --threads $N_CPU --force --output"
TABIX="tabix --threads $N_CPU"

# mask final junctions file 1
echo "masking junctions file 1"
zcat $SV_FINAL_JUNCTIONS_FILE_1 | 
perl ${MODULES_DIR}/analyze/SVs/package/package_jxns.pl |
$BGZIP $SV_MASKED_JUNCTIONS_FILE_1
checkPipe
$TABIX --sequence 1 --begin 2 --end 2 $SV_MASKED_JUNCTIONS_FILE_1
checkPipe

# mask final junctions file 2
echo "masking junctions file 2"
zcat $SV_FINAL_JUNCTIONS_FILE_2 | 
perl ${MODULES_DIR}/analyze/SVs/package/package_jxns.pl |
$BGZIP $SV_MASKED_JUNCTIONS_FILE_2
checkPipe
$TABIX --sequence 4 --begin 5 --end 5 $SV_MASKED_JUNCTIONS_FILE_2
checkPipe

echo
echo "done"

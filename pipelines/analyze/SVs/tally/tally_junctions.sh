# actions:
#   report some useful summary information to answer questions like:
#       what is the rate of SV junction artifacts per base?
# input:
#   $SV_ALIGNMENTS_FILE
#   $SV_FINAL_JUNCTIONS_FILE_1
# output:
#   various counts in a log report

echo "tallying SV junctions"

perl ${ACTION_DIR}/tally/tally_junctions.pl |
tee ${JUNCTION_TALLY_FILE}
checkPipe

echo "done"

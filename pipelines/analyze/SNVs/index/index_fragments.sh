# actions:
#   describe the SNV/indel combinations of single-alignment, i.e., non-SV, end-to-end RE fragments
#       in HiFiRe3, these are always PacBio HiFi reads, which are implicitly duplex and thus not stranded
#   use a streamlined, reduced-representation fragment index to count unique combinations
#   see index_fragments.pl for details
# input:
#   $SITE_SAM_PREFIX.$chrom.site_sam.gz, created by `analyze fragments`
# output:
#   $SNV_ALNS_PREFIX.<zeroPaddedChromIndex1>.txt.gz
#   $SNV_BASE_QUALITIES_FILE

# file paths
mkdir -p ${SNV_ALNS_DIR}
rm -f $SNV_ALNS_PREFIX.*.txt.gz

# index and parse unique fragment alignments
perl ${ACTION_DIR}/index/index_fragments.pl
checkPipe

# aggregate base quality information over all child threads
echo "aggregating base qualities"
cat <(
    echo -e "avgBaseQual\t*"
) <(
    cat $SNV_ALNS_PREFIX.base_qualities.*.txt | 
    sort -k1,1n | 
    bedtools groupby -g 1 -c 2 -o sum
) > $SNV_BASE_QUALITIES_FILE
checkPipe
cat $SNV_BASE_QUALITIES_FILE

echo
echo "done"

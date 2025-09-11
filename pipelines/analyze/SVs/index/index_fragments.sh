# actions:
#   describe the SV paths of all reads
#   use a streamlined, reduced-representation fragment index to count unique paths
#   see index_fragments.pl for details
# input:
#   $SITE_SAM_PREFIX.$chrom.site_sam.gz, created by `analyze fragments`
# output:
#   $SV_ALIGNMENTS_FILE
#   $SV_UNIQUE_JUNCTIONS_FILE
#   $SV_JUNCTION_SOURCES_FILE
#   library summation values (nReads, nBases...) to log file

# working variables
export SAMPLE_PATHS_PREFIX_WRK=${TMP_FILE_PREFIX_SMALL}.paths
rm -f $SAMPLE_PATHS_PREFIX_WRK.*.txt.bgz
rm -f $SAMPLE_PATHS_PREFIX_WRK.*.txt.gz
export SORT="sort --parallel=$N_CPU --buffer-size=4G --temporary-directory $TMP_DIR_WRK_SMALL"
PIGZ="pigz --stdout --processes $N_CPU"
BGZIP="bgzip --threads $N_CPU --force --output"
TABIX="tabix --threads $N_CPU --sequence 1"

# index and parse unique fragment alignments and junctions
if [[ "$CHECK_ENDPOINT_RE_MATCH" != "" && "$ENZYME_NAME" != "NA" ]]; then
    echo "indexing alignments by $ENZYME_NAME outer RE sites, i.e., ligFree"
    perl ${ACTION_DIR}/index/index_fragments.pl
    checkPipe
    SUMMATION_SCRIPT=${ACTION_DIR}/index/sum_not_deduplicated.sh
    if [ "$SEQUENCING_PLATFORM" == "ONT" ]; then
        SUMMATION_SCRIPT=${ACTION_DIR}/index/sum_deduplicated.sh
    fi
else 
    echo "indexing alignments by binned outer endpoints, i.e., tagFree or other" # may still perform junction site matching
    perl ${ACTION_DIR}/index/non_RE/index_fragments.pl  # the non-RE index works by chromosome, not site
    checkPipe
    SUMMATION_SCRIPT=${ACTION_DIR}/index/sum_deduplicated.sh
fi

# concenate pseudo-reference sequences into a single file for use in app
echo "concatenating read paths"
zcat $SAMPLE_PATHS_PREFIX_WRK.read_paths.*.txt.gz | 
$SORT -k1,1 -k3,3n | 
$BGZIP $SV_READ_PATHS_FILE
checkPipe
echo "indexing read paths" # indexed by QNAME and channel as a convenient way of indexed read retrieval
$TABIX --begin 3 --end 3 --zero-based $SV_READ_PATHS_FILE
checkPipe

# summarize the library read and base counts, with deduplication as needed
echo "calculating (deduplicated) library read and base counts"
zcat $SAMPLE_PATHS_PREFIX_WRK.alignment_sizes.*.txt.gz | 
bash $SUMMATION_SCRIPT | 
perl ${ACTION_DIR}/index/summarize_library.pl

# finish the sort on distal alignments and merge into first alignments, which were sorted during path indexing
#   when done, all alignment segments in kept paths are present
#       one row per alignment segment
#       sorted by chromIndex1,refPos1_5
#   unique path identifiers allow downstream intersection to other alignments or junctions in the same paths
#   for SVs, unique paths do NOT account for SNPs/indels in the alignments, only their reference spans and orientations
#   counts indicate how many times each unique set of alignments+junctions, i.e., each unique path, was observed
#   counts are NOT deduplicated, so ONT/tagFree count may be inflated here (unlike base summation above or deduplicated junctions)
#   metadata fields support downstream filters and plot devices
#   bgzip+tabix allows rapid retrieval for plotting in genome browser
echo "sorting distal alignments into first alignments"
$SORT --merge -k1,1n -k2,2n <(
    zcat $SAMPLE_PATHS_PREFIX_WRK.first_alignments.*.txt.gz
) <(
    zcat $SAMPLE_PATHS_PREFIX_WRK.distal_alignments.*.txt.gz | 
    $SORT -k1,1n -k2,2n -k3,3n
) | 
$BGZIP $SV_ALIGNMENTS_FILE
checkPipe
echo "indexing sorted alignments"
$TABIX --begin 2 --end 2 $SV_ALIGNMENTS_FILE # NOT --begin 2 --end 3 because column 3 can be less than column 2 
checkPipe

# sort the unique instances of junctions
#   when done, all unique junctions in kept paths are present
#       one row per ordered junction as node1,node2,jxnType = the junction key WITHOUT junction properties
#       sorted by chromIndex1_1,refPos1_1,strandIndex0_1,chromIndex1_2,refPos1_2,strandIndex0_2
#   summed counts indicate how many times each junction key was observed across all matching paths
#   alignmentOffset+jxnBases = junction properties taken from the most frequently observed path that contained each junction key
#   unique path identifiers allow downstream intersection to other alignments or junctions in the same paths
#   metadata fields support downstream filters and plot devices
#   unlike alignments, this is not the final output for junctions, which still require
#       merging event metadata from junction_sources into junctions
#       fuzzy matching of junction nodes to further aggregate inexact junction matches
#       ONT duplex purging when junctions appeared twice on different strands in the same channel
#       PCR duplicate purging when junctions appeared on sheared (not RE-cleaved) molecules with the same outer endpoints
#       creation of a second reverse-sorted junction file for indexed retrieval of both junction nodes for independent local plotting
echo "sorting unique junctions"
zcat $SAMPLE_PATHS_PREFIX_WRK.unique_junctions.*.txt.gz | 
$SORT -k1,1n -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n -k7,7n -k8,8nr |
bedtools groupby -g 1,2,3,4,5,6,7 -c 8,9,10,11,12,13 -o sum,first,first,distinct,max,min |
$PIGZ > $SV_UNIQUE_JUNCTIONS_FILE
checkPipe

# sort and aggregate the individual events that gave rise to the unique junctions
#   when done, file has one line per unique ordered junction WITH junction properties:
#       sorted and keyed by node1,node2,jxnType,alnOffset,jxnBases
#       sequencing event metadata collapsed to lists or taken as best-worst values
#   each line carries information to:
#       select the aggregated events that match the best unique junction instance for merging into unique junctions
#       assess ONT channel duplexes across all sequencing events matching a junction key
#       report on the identity, sequence and quality of the original sequencing events
echo "sorting and grouping junction sources"
zcat $SAMPLE_PATHS_PREFIX_WRK.junction_sources.*.txt.gz | 
$SORT -k1,1n -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n -k7,7n -k8,8 -k9,9 |
bedtools groupby -g 1,2,3,4,5,6,7,8,9 \
                 -c 10,11,12,13,14,15,16,17,18,19 \
                 -o collapse,collapse,collapse,distinct,max,min,max,collapse,collapse,collapse | 
$PIGZ > $SV_JUNCTION_SOURCES_FILE
checkPipe

echo
echo "done"

# actions:
#   merge event metadata from junction_sources into junctions
#   perform fuzzy matching of junction nodes to each other to aggregate inexact junction matches
#   purge ONT duplexes when junctions appeared on different strands in the same channel
#   record number of supporting reads per junction and whether any/all were inside/outside the allowed insert size range
#   TODO: analyze junctions for content flags, e.g. repetitive elements?
#   create a second reverse-sorted junction file for indexed retrieval of both junction nodes for independent local plotting
# input:
#   $SV_UNIQUE_JUNCTIONS_FILE
#   $SV_JUNCTION_SOURCES_FILE
# output:
#   $SV_FINAL_JUNCTIONS_FILE_1
#   $SV_FINAL_JUNCTIONS_FILE_2

# working variables
export SORT_COMMAND="sort --parallel=$N_CPU --buffer-size=4G --temporary-directory $TMP_DIR_WRK_SMALL"
BGZIP="bgzip --threads $N_CPU --force --output"
TABIX="tabix --threads $N_CPU"

# prepare for duplicate purging as needed, otherwise just cut the fields used for that purpose
PURGE_DUPLICATES="cut -f 1-26"
if [ "$SEQUENCING_PLATFORM" == "ONT" ]; then
    PURGE_DUPLICATES="perl ${ACTION_DIR}/aggregate/purge_ONT_duplexes.pl"
elif [ "$DEDUPLICATE_READS" == "TRUE" ]; then
    PURGE_DUPLICATES="perl ${ACTION_DIR}/aggregate/purge_PCR_duplicates.pl"
fi

# merge event metadata from junction_sources into junctions
perl ${ACTION_DIR}/aggregate/merge_events.pl |

# sort the junctions by chromIndex1_1,strandIndex0_1,refPos1_1 in preparation for breakpoint node1 fuzzy matching
$SORT_COMMAND -k1,1n -k3,3n -k2,2n |

# perform fuzzy matching of junction nodes to each other to aggregate inexact junction matches in SV groups
# when done, each row contains one fully filtered and aggregated unique output junction
perl ${ACTION_DIR}/aggregate/group_nodes.pl |

# purge ONT duplexes when junctions appeared on different strands in the same channel
# this action modifies junction counts in place in streamed junction rows
$PURGE_DUPLICATES |

# TODO: here, either within cross_reference.pl or as a separate step, 
# cross-compare called junctions to SVs expected by a sample's known genotype
# for now, just initialize the expected flag to 0
awk 'BEGIN{ OFS="\t" }{ print $0, 0; }' |

# re-align junctions to each other to try to resolve some artifact singletons
perl ${ACTION_DIR}/aggregate/cross_reference.pl |

# add final junction metadata (SV size, etc.)
# report junction coverage distribution to log stream
# repeat junctions in format suitable for breakpoint coverage assessment
perl ${ACTION_DIR}/aggregate/junction_coverage.pl |

# add coverage information for each junction breakpoint (all reads, not just SV reads)
# repeat junctions for saving to disk
bedtools intersect -loj -a <(
    zcat $SV_ALIGNMENTS_FILE | 
    awk 'BEGIN{OFS="\t"}{ 
        if ($2 <= $3) {
            print $1, $2 - 1, $3, $10;
        } else {
            print $1, $3 - 1, $2, $10;
        }
    }'
) -b stdin |
awk '$6 >= 0' | # since intersect --loj reports null B features with -1 positions
perl ${ACTION_DIR}/aggregate/breakpoint_coverage.pl |

# sort and index the junctions by node1
$SORT_COMMAND -k1,1n -k2,2n -k3,3n -k4,4n -k5,5n -k6,6n|
$BGZIP $SV_FINAL_JUNCTIONS_FILE_1
checkPipe
$TABIX --sequence 1 --begin 2 --end 2 $SV_FINAL_JUNCTIONS_FILE_1
checkPipe

# resort and index the junctions file by node2
zcat $SV_FINAL_JUNCTIONS_FILE_1 | 
$SORT_COMMAND -k4,4n -k5,5n -k6,6n -k1,1n -k2,2n -k3,3n | 
$BGZIP $SV_FINAL_JUNCTIONS_FILE_2
checkPipe
$TABIX --sequence 4 --begin 5 --end 5 $SV_FINAL_JUNCTIONS_FILE_2
checkPipe

echo
echo "done"

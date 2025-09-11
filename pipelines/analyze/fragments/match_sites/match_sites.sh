# actions:
#   master script that prepares aligned fragments for downstream variant analysis
#     purge unmapped reads
#     split anomalous sequence gaps into paired orphaned reads
#     order alignments across each read from 5' to 3'
#     renumber read2 orphans to read1
#     repack SAM to SITE_SAM for streamlined parsing
#     calculate the traversal delta over all sets of read alignments (across all sets of junctions) to set blockN
#         used to suppress spurious junctions arising in low quality regions with ~equal inserted and deleted bases
#     enforce alignment-level quality filtering with read truncation/cutting in 5' to 3' order
#         suppresses untrusted junctions arising from aligner issues
#     match alignment nodes to RE filtering sites (tagFree and ligFree)
#     discard sequences where one or both complete outer nodes fail to match an RE filtering site (ligFree only)
#     split chimeras at SV breakpoint nodes that match (discard distal/3' portion):
#         - RE sites (due to intermolecular ligation; tagFree and ligFree)
#         - foldback inversions (due to ONT duplex reads or inversion synthesis)
#         - low quality insertions (ONT only, to suppress follow-on chimeras)
#         - adapter insertions (due to any process that conjoins reads across an adapter)
#     project the 3' end of required alignments to the next (haplotype-consistent) RE site (ligFree only)
#   no site matching is performed when analyzing libraries made from random, not RE, sheared gDNA
# input:
#   NAME_BAM_FILE, must have SEQ/QUAL set to * when not needed, with xf:i: and xh:i tags
#   RE_SITES_DIR, as populated by `genome index` or `locate` (tagFree and ligFree, NA otherwise)
#   library-level config parameters:
#       CHECK_ENDPOINT_RE_MATCH (true for ligFree only, where endpoints should match RE sites)
#       CHECK_JUNCTION_RE_MATCH (true for tagFree and ligFree, where junction nodes should not match RE sites)
#   platform-level config parameters that determine matching tolerance
#       ACCEPT_ENDPOINT_DISTANCE, stricter, for ligFree outer endpoints
#       REJECT_JUNCTION_DISTANCE, more permissive, thus rejects junction chimera more easily
# output:
#   SITE_SAM written to chromosome-level files in SITE_SAM_DIR

# variables
MATCH_DIR=${ACTION_DIR}/match_sites
UNMAPPED=4
if [[ "$CHECK_JUNCTION_RE_MATCH" != "" && "$ENZYME_NAME" != "NA" ]]; then # always true if CHECK_ENDPOINT_RE_MATCH is true
    MATCH_NODES_SCRIPT=${MATCH_DIR}/match_nodes_to_sites.pl # outer nodes are matched, somewhat pointlessly, for tagFree
    SPLIT_CHIMERAS_COMMAND="perl ${MATCH_DIR}/split_chimeric_reads.pl"
    LOAD_RE_LOOKUPS=TRUE
else 
    MATCH_NODES_SCRIPT=${MATCH_DIR}/non_RE/fill_site_pos.pl # but "other" non-RE libraries don't match at all
    SPLIT_CHIMERAS_COMMAND="cat" 
fi
if [[ "$CHECK_ENDPOINT_RE_MATCH" != "" && "$ENZYME_NAME" != "NA" ]]; then
    DISCARD_UNMATCHED_SCRIPT=${MATCH_DIR}/discard_unmatched_sequences.pl # outer endpoint matching and projection for ligFree only
    PROJECT_READS_SCRIPT=${MATCH_DIR}/project_incomplete_reads.pl
else 
    DISCARD_UNMATCHED_SCRIPT=${MATCH_DIR}/non_RE/set_end_to_end.pl
    PROJECT_READS_SCRIPT=${MATCH_DIR}/non_RE/finish_non_RE_reads.pl
fi

# load the site lookup binaries into /dev/shm for fastest scanning
# TODO: load several sites to support admixed ligFree libraries? (requires confident assignment of RE source of all fragments)
# production
if [ "$LOAD_RE_LOOKUPS" == "TRUE" ]; then
    echo "copying $ENZYME_NAME site lookup files to ram disk"
    cp $CLOSEST_SITE_LOOKUP    $TMP_DIR_WRK_SHM
    cp $SITE_DATA_LOOKUP   $TMP_DIR_WRK_SHM
    export CLOSEST_SITE_LOOKUP_WRK=$TMP_DIR_WRK_SHM/`basename $CLOSEST_SITE_LOOKUP`
    export SITE_DATA_LOOKUP_WRK=$TMP_DIR_WRK_SHM/`basename $SITE_DATA_LOOKUP`
    # debugging
    # export CLOSEST_SITE_LOOKUP_WRK=$CLOSEST_SITE_LOOKUP
    # export SITE_DATA_LOOKUP_WRK=$SITE_DATA_LOOKUP
fi
echo "parsing alignments to quality-filtered, RE-matched, error-corrected fragments and reads"

# working per alignment:
#   purge unmapped reads/empty alignments and convert to SAM
#   orphaned paired reads continue on, possibly as isolated read2 for now
#   strip to just the tags needed for fragment parsing, which must be in order de,cs,xf,xh
samtools view --threads $SAMTOOLS_CPU --exclude-flags $UNMAPPED --keep-tag de,cs,xf,xh $NAME_BAM_FILE | 

# working per event:
#   PENDING: add QNAME extensions to external alignments
#   enforce zoo rejection
#   split anomalous gaps: append splitGapReadN to all alignment QNAME
#   renumber read2 orphans as read1
#   order alignments across each read from 5' to 3'
#   repack SAM to SITE_SAM for streamlined parsing
#   record number of aligned read bases and associated ref span
perl ${MATCH_DIR}/order_alignments.pl |

# working per read:
#   calculate the reference vs. read traversal delta between all sets of alignments (across all sets of junctions)
#   use calculated deltas and MIN_TRAVERSAL_DELTA to set BLOCK_N of all alignments
#   junctions between alignments with the same BLOCK_N may be considered artifactual during SV analysis
#   block numbering is done before alignment quality truncation when all read alignments are still present
#   SVs that fail block traversal checks may or may not be split by downstream quality filters
perl ${MATCH_DIR}/check_traversal_delta.pl | 

# working per sequence:
#   check whether alignments pass alignment-level quality filters
#   when an alignment fails, remove it and all alignments further 3' on the read
#   renumber new read2 orphans as read1
#   must come after order_alignments as quality filtering is executed from 5' to 3' on read
perl ${MATCH_DIR}/enforce_quality_filters.pl |

# working per alignment:
#   find the nearest filtering RE site to each alignment node and the direction to it
#   this analysis is haplotype-independent (but looks up site haplotypes)
perl ${MATCH_NODES_SCRIPT} |

# working per sequence:
#   reject and discard sequences that don't match RE sites at complete outermost endpoints
#       disregard incomplete 3' ends when applying this filter
#       this filter is haplotype-independent, i.e., any site on any haplotype is accepted
#   fill the sequence 3' site fields end when already known for end-to-end/complete sequences with no SV
#   append readN to QNAME to force all subsequent actions to act per read
perl ${DISCARD_UNMATCHED_SCRIPT} | 

# working per read:
#   split SV junctions that match RE sites or are otherwise identified as chimeras
#   in particular, reads are also processed for adapter presence, ONT follow-ons, and foldback inversion chimeras
#   discard the distal/3' portion of split chimeric reads as unusable, retaining the 5' end
${SPLIT_CHIMERAS_COMMAND} | 

# working per read:
#   when not set above, project the 3' end of required alignments to the next (haplotype-consistent) RE site
#   assess the proper fragment insert size distribution
#   assess whether alignment/read ends matched a target region, if TARGETS_BED provided
#   print SITE_SAM to chromosome-level file in preparation for variant analysis
perl ${PROJECT_READS_SCRIPT}

# action:
#     set environment variables for performing variant analysis
# requires:
#     source $MODULES_DIR/library/set_library_vars.sh
# usage:
#     source $MODULES_DIR/analyze/set_analysis_vars.sh

# analysis identifiers
export ANALYSIS_IDENTIFIER=${DATA_NAME}.analysis
export ANALYSIS_PREFIX=${TASK_DIR}/${ANALYSIS_IDENTIFIER}

# fragment analysis output
export SITE_SAM_DIR=${TASK_DIR}/site_sam
mkdir -p ${SITE_SAM_DIR}
export SITE_SAM_PREFIX=${SITE_SAM_DIR}/${ANALYSIS_IDENTIFIER}
export ANALYSIS_CHROMS_FILE=${ANALYSIS_PREFIX}.chroms.txt

# post-filtering insert sizes file, as created by analyze/fragments, with filterd and projected sizes
export FILTERED_INSERT_SIZES_FILE=${DATA_FILE_PREFIX}.insert_sizes.filtered.txt
export FILTERED_STEM_LENGTHS_FILE=${DATA_FILE_PREFIX}.stem_lengths.filtered.txt

# SV analysis output
export SV_READ_PATHS_FILE=${ANALYSIS_PREFIX}.read_paths.txt.bgz # one line per read with a map of all junctions
export SV_ALIGNMENTS_FILE=${ANALYSIS_PREFIX}.alignments.txt.bgz
export SV_UNIQUE_JUNCTIONS_FILE=${ANALYSIS_PREFIX}.unique_junctions.txt.gz
export SV_JUNCTION_SOURCES_FILE=${ANALYSIS_PREFIX}.junction_sources.txt.gz
export SV_FINAL_JUNCTIONS_FILE_1=${ANALYSIS_PREFIX}.final_junctions_1.txt.bgz # sorted and indexed on node1
export SV_FINAL_JUNCTIONS_FILE_2=${ANALYSIS_PREFIX}.final_junctions_2.txt.bgz # redundant with SV_FINAL_JUNCTIONS_FILE_1, sorted and indexed on node2
export SV_MASKED_JUNCTIONS_FILE_1=${ANALYSIS_PREFIX}.final_junctions_1.masked.txt.bgz # as above, but with SEQ, QUAL and CIGAR masked to *
export SV_MASKED_JUNCTIONS_FILE_2=${ANALYSIS_PREFIX}.final_junctions_2.masked.txt.bgz 
export JUNCTION_TALLY_FILE=${DATA_GENOME_PREFIX}.junction_tally.txt
export GENOME_TALLY_FILE=${DATA_GENOME_PREFIX}.genome_tally.txt

# SNV analysis output
export SNV_ALNS_DIR=${TASK_DIR}/snv_alns
export SNV_ALNS_PREFIX=${SNV_ALNS_DIR}/${DATA_NAME}.${GENOME}.snv_alns
# export SNV_ALIGNMENT_MAPS_FILE=${ANALYSIS_PREFIX}.alignment_maps.txt.bgz
export SNV_BASE_QUALITIES_FILE=${ANALYSIS_PREFIX}.base_qualities.txt
export SNV_CHROM_FILE_PREFIX=${SNV_ALNS_DIR}/${DATA_NAME}.${GENOME} # initial chrom-level files (gzipped when applicable)
export SNV_GENOME_PILEUP_PREFIX=${SNV_CHROM_FILE_PREFIX}.genome_pileup
export SNV_READ_PILEUP_PREFIX=${SNV_CHROM_FILE_PREFIX}.read_pileup
export SNV_SUMMARY_TABLE_PREFIX=${SNV_CHROM_FILE_PREFIX}.snv_summary
export SNV_GENOME_PILEUP=${DATA_GENOME_PREFIX}.genome_pileup.bed.bgz # final aggregated genome/read-level files (bgzipped and tabixed when applicable)
export SNV_READ_PILEUP=${DATA_GENOME_PREFIX}.read_pileup.txt.gz
export SNV_SUMMARY_TABLE=${DATA_GENOME_PREFIX}.snv_summary.bed.bgz

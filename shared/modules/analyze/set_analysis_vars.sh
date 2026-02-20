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
export ANALYSIS_CHROMS_FILE=${ANALYSIS_PREFIX}.chroms.txt

# post-filtering insert sizes file, as created by analyze/fragments, with filterd and projected sizes
export FILTERED_INSERT_SIZES_FILE=${DATA_FILE_PREFIX}.insert_sizes.filtered.txt
export FILTERED_STEM_LENGTHS_FILE=${DATA_FILE_PREFIX}.stem_lengths.filtered.txt

# SV analysis output
export SV_SAMPLES_FILE=${ANALYSIS_PREFIX}.sv_samples.txt
export SV_READ_PATHS_FILE=${ANALYSIS_PREFIX}.read_paths.txt.bgz # one line per read with a map of all junctions
export SV_ALIGNMENTS_FILE=${ANALYSIS_PREFIX}.alignments.txt.bgz
export SV_COVERAGE_FILE=${ANALYSIS_PREFIX}.coverage.bed.bgz
export SV_FINAL_JUNCTIONS_FILE_1=${ANALYSIS_PREFIX}.final_junctions_1.txt.bgz # sorted and indexed on node1
export SV_FINAL_JUNCTIONS_FILE_2=${ANALYSIS_PREFIX}.final_junctions_2.txt.bgz # redundant with SV_FINAL_JUNCTIONS_FILE_1, sorted and indexed on node2
export JUNCTION_TALLY_FILE=${DATA_GENOME_PREFIX}.junction_tally.txt
export GENOME_TALLY_FILE=${DATA_GENOME_PREFIX}.genome_tally.txt

# SNV analysis output
export SNV_SAMPLES_FILE=${ANALYSIS_PREFIX}.snv_samples.txt
export SNV_ALL_READS_PILEUP_BGZ=${DATA_GENOME_PREFIX}.all_reads.pileup.bed.bgz 
export SNV_ALL_READS_VARIANTS_BGZ=${DATA_GENOME_PREFIX}.all_reads.snv_indel.txt.bgz
export SNV_ERROR_CORRECTED_PILEUP_BGZ=${DATA_GENOME_PREFIX}.error_corrected_${MIN_N_PASSES}.pileup.bed.bgz 
export SNV_ERROR_CORRECTED_VARIANTS_BGZ=${DATA_GENOME_PREFIX}.error_corrected_${MIN_N_PASSES}.snv_indel.txt.bgz

# action:
#     set environment variables to guide a subsequent read sequence alignment
# expects:
#     source $MODULES_DIR/genome/set_genome_vars.sh
# usage:
#     source $MODULES_DIR/align/set_alignment_vars.sh

# set the alignment log files
export FASTP_LOG_PREFIX=$LOG_FILE_PREFIX.fastp
export MINIMAP_LOG_FILE=$LOG_FILE_PREFIX.minimap.log

# set the product alignment files
#   HTSlib has a limit of 268,435,456 bases for any single CIGAR operation (soft clips are the main concern)
#   this creates errors in BAM format, so must enforce read length limit
#   this is not usually a concern since RE fragments are almost never that long
#   CRAM does not have this concern, but name-sorted CRAM files are inefficient and discouraged
export NAME_BAM_FILE=$DATA_GENOME_PREFIX.name.bam
export COORDINATE_BAM_FILE=$DATA_GENOME_PREFIX.coordinate.bam
export COORDINATE_BAM_INDEX=$COORDINATE_BAM_FILE.bai
export NAME_REALIGNED_BAM_FILE=$DATA_GENOME_PREFIX.name.realigned.bam
export COORDINATE_REALIGNED_BAM_FILE=$DATA_GENOME_PREFIX.coordinate.realigned.bam

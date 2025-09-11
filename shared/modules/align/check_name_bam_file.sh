# action:
#     parse and check the input name-sorted bam file
# expects:
#     source $MODULES_DIR/align/set_alignment_vars.sh
# optional:
#     user-provided $NAME_BAM, to override default $NAME_BAM_FILE

# allow users to provide their own name-sorted bam file
export IS_USER_BAM=0
if [[ "$NAME_BAM" != "" && "$NAME_BAM" != "null" ]]; then
    export IS_USER_BAM=1
    export NAME_BAM_FILE=$NAME_BAM # can also be a BAM file
fi

# make sure the input bam file exists
if [ ! -e $NAME_BAM_FILE ]; then
    echo -e "bam file not found:\n    $NAME_BAM_FILE"
    exit 1
fi

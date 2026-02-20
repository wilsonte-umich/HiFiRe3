# action:
#     set environment variables for multi-sample comparisons
# requires:
#     source $MODULES_DIR/library/set_library_vars.sh
# usage:
#     source $MODULES_DIR/compare/set_compare_vars.sh

export NAME_BAM_FILES=""
if [ "$COMPARE_SAMPLE_DIRS" != "null" ]; then
    COMPARE_SAMPLE_DIRS=$(echo $COMPARE_SAMPLE_DIRS | tr ',' ' ')
else
    if [ "$COMPARE_PROJECT_DIR" == "null" ]; then
        export COMPARE_PROJECT_DIR=$TASK_DIR
    fi
    COMPARE_SAMPLE_DIRS=$(ls -d $COMPARE_PROJECT_DIR/*/)
fi

export NAME_BAM_FILES=""
for COMPARE_SAMPLE_DIR in $COMPARE_SAMPLE_DIRS; do
    SAMPLE_NAME_BAM_FILE=$COMPARE_SAMPLE_DIR/*.name.bam
    if [ -f $SAMPLE_NAME_BAM_FILE ]; then
        SAMPLE_NAME_BAM_FILE=$(ls $SAMPLE_NAME_BAM_FILE)
        if [ "$NAME_BAM_FILES" == "" ]; then
            export NAME_BAM_FILES="$SAMPLE_NAME_BAM_FILE"
        else
            export NAME_BAM_FILES="$NAME_BAM_FILES,$SAMPLE_NAME_BAM_FILE"
        fi
    fi
done

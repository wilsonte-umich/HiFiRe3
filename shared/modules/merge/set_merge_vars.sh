# action:
#     set environment variables for multi-sample junction list merges
# usage:
#     source $MODULES_DIR/merge/set_merge_vars.sh

if [ "$MERGE_SAMPLE_DIRS" != "null" ]; then
    MERGE_SAMPLE_DIRS=$(echo $MERGE_SAMPLE_DIRS | tr ',' ' ')
else
    if [ "$MERGE_PROJECT_DIR" == "null" ]; then
        export MERGE_PROJECT_DIR=$TASK_DIR
    fi
    MERGE_SAMPLE_DIRS=$(ls -d $MERGE_PROJECT_DIR/*/)
fi

export MERGE_INPUT_DIRS=""
for MERGE_SAMPLE_DIR in $MERGE_SAMPLE_DIRS; do
    SAMPLE_FINAL_JUNCTION_FILE=$MERGE_SAMPLE_DIR/*.final_junctions_1.txt.bgz
    if [ -f $SAMPLE_FINAL_JUNCTION_FILE ]; then
        if [ "$MERGE_INPUT_DIRS" == "" ]; then
            export MERGE_INPUT_DIRS="$MERGE_SAMPLE_DIR"
        else
            export MERGE_INPUT_DIRS="$MERGE_INPUT_DIRS,$MERGE_SAMPLE_DIR"
        fi
    fi
done

N_MERGE_INPUT_DIRS=$(echo $MERGE_INPUT_DIRS | tr ',' ' ' | wc -w)
if [ "$N_MERGE_INPUT_DIRS" == "0" ]; then
    echo "no final junction files found to merge"
    exit 1
fi
if [ "$N_MERGE_INPUT_DIRS" == "1" ]; then
    echo "only one final junction file found, nothing to merge"
    exit 1
fi

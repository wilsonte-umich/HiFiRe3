# developer script to create options documentation files
# run from the folder containing this script

TOOL_SUITE="HiFiRe3"
SUITE_COMMAND="hf3"

MDI="../../../../mdi"

create_options_file () {

    local FILE="${PIPELINE}_${ACTION}.md"
    local COMMAND="$PIPELINE $ACTION --help"

    echo $PIPELINE $ACTION
    
    echo -e "## $PIPELINE $ACTION options\n" > $FILE

    echo -e "|Tool Suite|Pipeline|Action|" >> $FILE
    echo -e "|---|---|---|" >> $FILE
    echo -e "|$TOOL_SUITE|$PIPELINE|$ACTION|\n" >> $FILE

    echo -e '```'"\n$SUITE_COMMAND $COMMAND\n"'```'"\n" >> $FILE

    echo -e '```' >> $FILE
    $MDI -d $COMMAND | tail -n+2 >> $FILE
    echo -e '```' >> $FILE
}

PIPELINES=`ls -1 ../pipelines | grep -v _`
for PIPELINE in $PIPELINES; do
    PIPELINE_DIR="../pipelines/$PIPELINE"
    if [ -d $PIPELINE_DIR ]; then
        ACTIONS=`ls -1 $PIPELINE_DIR | grep -v _`
        for ACTION in $ACTIONS; do
            ACTION_DIR="$PIPELINE_DIR/$ACTION"
            if [ -d $ACTION_DIR ] && [ -f $ACTION_DIR/Workflow.sh ]; then
                create_options_file
            fi
        done
    fi
done

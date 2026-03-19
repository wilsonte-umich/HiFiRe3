#!/bin/bash

# set read file variables and paths
source $MODULES_DIR/rust/set_rust_vars.sh
source $MODULES_DIR/align/set_read_file_vars.sh

# create backup directory for original files
# can be deleted en bloc later after confirming reformatting
BU_DIR=${READ_FILE_DIR}/bu_original_files
mkdir -p ${BU_DIR}

echo "reformatting unaligned BAM files in directory"
echo "  $READ_FILE_DIR"
for UBAM_FILE in $READ_1_FILES; do
    FILE_NAME=$(basename $UBAM_FILE)
    echo "    $FILE_NAME"

    # move original file to backup directory
    echo "      moving to backup directory"
    BU_FILE=${BU_DIR}/${FILE_NAME}
    mv ${UBAM_FILE} ${BU_FILE}

    # run the reformatting
    echo "      reformatting and compressing"
    samtools view -@ ${N_CPU} --with-header ${BU_FILE} |
    ${HF3_TOOLS_BIN} reformat_ont |
    samtools view -@ ${N_CPU} -b -o ${UBAM_FILE} -
    checkPipe

    # delete the backup (pending)
    echo "      done (backup file NOT deleted, see .../bu_original_files)"

done

echo "reformatting complete"

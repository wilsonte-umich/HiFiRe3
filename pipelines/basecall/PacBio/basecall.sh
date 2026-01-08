#!/bin/bash

# run basecalling; all file handling is autonomous within the tool
if [[ "${FORCE_BASECALLING}" != "" && "${FORCE_BASECALLING}" != "0" ]]; then FORCE_BASECALLING="true"; fi
if [[ -f "${READ_BAM_FILE}" && "${FORCE_BASECALLING}" != "true" ]]; then
    echo "Basecalled read BAM file ${READ_BAM_FILE} already exists; skipping basecalling step."
    echo "Set flag --force-basecalling to re-run basecalling."
    exit 0
fi
$SUITE_BIN_DIR/hf3_tools basecall_pacbio
checkPipe

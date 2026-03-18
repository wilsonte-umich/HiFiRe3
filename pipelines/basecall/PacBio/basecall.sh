#!/bin/bash

# run basecalling; all file handling is autonomous within the tool
if [[ "${FORCE_BASECALLING}" != "" && "${FORCE_BASECALLING}" != "0" ]]; then FORCE_BASECALLING="true"; fi
if [[ -f "${READ_BAM_FILE}" && "${FORCE_BASECALLING}" != "true" ]]; then
    echo "Basecalled read BAM file ${READ_BAM_FILE} already exists; skipping basecalling step."
    echo "Set flag --force-basecalling to re-run basecalling."
    exit 0
fi

# as needed, create and save the appropriate minimap2 genome alignment index
export GENOME_FASTA_WRK=${GENOME_FASTA}
export ALIGNMENT_MODE_WRK="map-hifi"
export CREATE_MM2_INDEX_MESSAGE="continuing with basecalling"
source ${MODULES_DIR}/align/create_mm2_index.sh # sets variable ${MINIMAP2_INDEX_WRK}

# run basecalling
${HF3_TOOLS_BIN} basecall_pacbio
checkPipe

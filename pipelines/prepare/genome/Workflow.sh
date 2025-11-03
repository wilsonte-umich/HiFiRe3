#!/bin/bash

# function to display alert messages
show_alert_message () {
    echo
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo -e "$1"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo
}

# check that the requested genome is supported
export SUPPORTED_GENOMES="hs1|CHM13|hg38|GRCh38|mm39|GRCm39|dm6"
if [[ "|${SUPPORTED_GENOMES}|" != *"|${GENOME}|"* ]]; then
    show_alert_message \
"genome not supported: $GENOME\n"\
"supported genomes names are: $SUPPORTED_GENOMES"
    exit 1
fi

# coerce GRC to UCSC genome names
# using UCSC names helps support the MDI genome browser
#   which dynamically reads genome metadata from the UCSC API
export UCSC_GENOMES="hs1|hg38|mm39|dm6"
if [[ "|${UCSC_GENOMES}|" != *"|${GENOME}|"* ]]; then
    GRC_GENOME=${GENOME}
    if [ "$GENOME" = "CHM13" ]; then
        export GENOME=hs1
    elif [ "$GENOME" = "GRCh38" ]; then
        export GENOME=hg38
    elif [ "$GENOME" = "GRCm39" ]; then
        export GENOME=mm39
    fi
    show_alert_message \
"GRC genome name ${GRC_GENOME} was coerced to the UCSC genome name ${GENOME}\n"\
"usage is optimized for UCSC genome names as they directly support the genome browser"
fi

# set derivative environment variables and file paths
export GENOME_DIR=${GENOMES_DIR}/${GENOME}
mkdir -p ${GENOME_DIR}
export HIFIRE3_PREPARING_GENOME=TRUE
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/REs/set_re_vars.sh

# download the genome fasta and supporting metadata
runWorkflowStep 1 download download.sh

# update derivative environment variables with genome metadata
export HIFIRE3_PREPARING_GENOME=FALSE
source ${MODULES_DIR}/genome/set_genome_vars.sh
source ${MODULES_DIR}/REs/set_re_vars.sh

# collect all genome RE sites for a set of restriction enzymes
runWorkflowStep 2 digest ${MODULES_DIR}/genome/digest.sh

# analyze the restriction fragment length distributions for plotting in app
runWorkflowStep 3 assemble ${MODULES_DIR}/genome/assemble.sh

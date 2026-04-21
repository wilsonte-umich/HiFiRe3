# action:
#     check for the existence of the required reference genome support files
#     if needed, perform the same actions as `prepare genome`
#     thus, inline genome preparation without requiring a separate submitted job
#     this script does not handle composite genomes, an advanced usage that requires running `combine genomes`
# expects:
#     $GENOME_DIR, will default to $MDI_DIR/resources/genomes/$GENOME
#     $GENOME
# usage:
#     in Workflow.sh:
#         runWorkflowStep 1 check_genome $MODULES_DIR/genome/check_genome.sh
#     prior to calling 
#         source $MODULES_DIR/genome/set_genome_vars.sh

# set variables to the required genome files
export HIFIRE3_PREPARING_GENOME=TRUE
source ${MODULES_DIR}/genome/parse_genome.sh
source ${MODULES_DIR}/genome/set_genome_vars.sh

# check genome file existence; download single genomes if missing
# some files may be empty but are always created by the downloader
GENOME_FASTA_LOCK=${GENOME_FASTA}.lock
while [ -f ${GENOME_FASTA_LOCK} ]; do # prevent two jobs from creating the same file at the same time
    echo "waiting for $GENOME genome to be downloaded"
    sleep 1m
done
if [[ ! -f ${GENOME_FASTA} || ! -f ${GENOME_GAPS_FILE} || ! -f ${GENOME_EXCLUSIONS_BED} || ! -f ${GENES_BED} ]]; then
    if [[ "${IS_COMPOSITE_GENOME}" != "" &&  "${IS_COMPOSITE_GENOME}" != "0" ]]; then
        echo "one or more missing composite files for genome $GENOME"
        echo "create them using the `combine genomes` pipeline action"
        exit 1
    else
        touch ${GENOME_FASTA_LOCK}
        source ${MODULES_DIR}/genome/download.sh
        rm ${GENOME_FASTA_LOCK}
    fi
fi

# if needed, check RE site file existence; perform in silico digestion if missing
if [[ "${ENZYME_NAME}" != "NA" && "${ENZYME_NAME}" != "null" && "${ENZYME_NAME}" != "" ]]; then
    GENOME_SITES_BGZ_LOCK=${GENOME_SITES_BGZ}.lock
    while [ -f ${GENOME_SITES_BGZ_LOCK} ]; do # prevent two jobs from creating the same file at the same time
        echo "waiting for $GENOME genome to be digested in silico"
        sleep 5m
    done
    export MIN_PRIORITY_LEVEL=4
    source ${MODULES_DIR}/REs/set_re_vars.sh
    if [ ! -f  ${GENOME_SITES_BGZ} ]; then
        touch ${GENOME_SITES_BGZ_LOCK}
        export HIFIRE3_PREPARING_GENOME=FALSE
        source ${MODULES_DIR}/genome/set_genome_vars.sh
        echo "performing in silico digestion of ${GENOME}"
        source ${MODULES_DIR}/genome/digest.sh
        source ${MODULES_DIR}/genome/assemble.sh
        rm ${GENOME_SITES_BGZ_LOCK}
    fi
fi

# reset for calling code
export HIFIRE3_PREPARING_GENOME=FALSE
echo "all required genome files are present"

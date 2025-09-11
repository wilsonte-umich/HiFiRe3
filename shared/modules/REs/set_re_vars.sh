# action:
#     set environment variables for known restriction enzymes
# requires:
#     source $MODULES_DIR/genome/set_genome_vars.sh
# usage:
#     source $MODULES_DIR/REs/set_re_vars.sh

# compatible restriction enyzmes
export BLUNT_RE_TABLE=${MODULES_DIR}/REs/blunt_enzymes.csv

# genome+RE-specific paths
export GENOME_REMAPS_DIR=${GENOME_DIR}/RE_maps
export GENOME_SITES_GZ=`echo ${GENOME_REMAPS_DIR}/${GENOME}.digest.${ENZYME_NAME}_*.txt.gz`
export GENOME_ENZYME_DIR=${GENOME_DIR}/${GENOME}_${ENZYME_NAME}
export GENOME_ENZYME_PREFIX=${GENOME_ENZYME_DIR}/${GENOME}.${ENZYME_NAME}


# set the dorado directories
echo "checking nanopore resource directories"
NANOPORE_DIR=${MDI_DIR}/resources/nanopore
DORADO_DIR=${NANOPORE_DIR}/dorado
mkdir -p ${DORADO_DIR}
DORADO_VERSION_NAME=dorado-${DORADO_VERSION}
DORADO_VERSION_DIR=${DORADO_DIR}/${DORADO_VERSION_NAME}
DORADO_EXECUTABLE=${DORADO_VERSION_DIR}/bin/dorado
DORADO_MODELS_DIR=${DORADO_DIR}/models
mkdir -p ${DORADO_MODELS_DIR}

# obtain dorado if not already present
if [ ! -d ${DORADO_VERSION_DIR} ]; then
    echo "downloading Dorado version "${DORADO_VERSION}
    TMP_DIR=${PWD}
    cd ${DORADO_DIR}
    DORADO_ARCHIVE=${DORADO_VERSION_NAME}".tar.gz"
    wget --no-verbose "https://cdn.oxfordnanoportal.com/software/analysis/"${DORADO_ARCHIVE}
    echo "extracting Dorado"
    tar -xzf ${DORADO_ARCHIVE}
    rm -f ${DORADO_ARCHIVE}
    cd ${TMP_DIR}
fi
if [ ! -f ${DORADO_EXECUTABLE} ]; then
    echo "missing Dorado executable: "${DORADO_EXECUTABLE}
    exit 1
fi 

# process and set input and output paths
echo "parsing input and output paths"
EXPANDED_INPUT_DIR=`echo ${POD5_DIR}` # could be a tar archive at this point
UBAM_DIR=${TASK_DIR}/ubam

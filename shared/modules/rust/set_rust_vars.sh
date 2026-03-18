# actions:
#   set the path to the hf3_tools binary based on the code version in use
#   download the binary from GitHub as needed
# input:
#   $SUITE_VERSION
#   $DEVELOPER_MODE
# output:
#   $HF3_TOOLS_BIN

getVersionedBinary wilsontelab/${SUITE_NAME} hf3_tools
export HF3_TOOLS_BIN=${VERSIONED_BINARY_PATH}

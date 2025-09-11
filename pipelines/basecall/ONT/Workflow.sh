#!/bin/bash

# set derivative environment variables and file paths
source ${MODULES_DIR}/REs/set_re_vars.sh
source $MODULES_DIR/agFree/set_library_vars.sh

# create temp directories
if [ "$POD5_BUFFER" = "shm" ]; then 
    source $MODULES_DIR/utilities/shell/create_temp_dir_shm.sh
    POD5_BUFFER_DIR=$TMP_DIR_WRK_SHM;   # prefer to use /dev/shm for files pod5/dorado actively use
else
    source $MODULES_DIR/utilities/shell/create_temp_dir_small.sh
    POD5_BUFFER_DIR=$TMP_DIR_WRK_SMALL; # but allow fall back to SSD if files too big for /dev/shm
fi

# parse MDI to Dorado options in preparation for basecalling
source $ACTION_DIR/parse_dorado_options.sh

# convert ONT POD5 event files to unaligned BAM, i.e., call bases
runWorkflowStep 1 basecall $ACTION_DIR/basecall_autobatch.sh

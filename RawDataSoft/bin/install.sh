#!/bin/bash

# where am I?
SCRIPT="$(readlink -f ${BASH_SOURCE[0]})"
DIRBIN="$(dirname $SCRIPT)"

# our software is here 
SOFTDIR="$(dirname $DIRBIN)"

# check if the environment for RAWDATASOFT exists already
if [ -z "$RAWDATASOFT" ]; then
    #echo "Setting up software environment for installation\n in $SOFTDIR"
    export PATH=$PATH:${SOFTDIR}/bin:${SOFTDIR}
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SOFTDIR}/lib
else
    #echo "Cleaning environment from existing installation ..."
    export PATH=$(echo $PATH | sed "s|$RAWDATASOFT|$SOFTDIR|g" | sed "s/::/:/g")
    export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | sed "s|$RAWDATASOFT|$SOFTDIR|g" | sed "s/::/:/g")
fi;

# soft directory
export RAWDATASOFT=${SOFTDIR}
# config directory

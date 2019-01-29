#!/bin/bash

ALGORITHM=$1
HEURISTIC=$2

MAX_TIME="30m"
ERRORS='errors.out'
if [ -z "$3" ]
then
    INSTANCE_PATH='./instances'
elif [ "$3" = "random" ]
then
    EXTRA=$3
else
    INSTANCE_PATH=$3
fi

INSTANCE_LIST=$(find ${INSTANCE_PATH} -type f -name '*.ideal')

function experiment {
    filename=$(basename -- $1)
    filename="${filename%.*}"
    filename="${filename}.test"
    timeout "${MAX_TIME}" sage experiment.sage $1 "${ALGORITHM}" "${HEURISTIC}" "${EXTRA}" > "${filename}" 2>> "${ERRORS}"
    if [ $? -eq 124 ]
    then
        echo "inf inf inf inf inf inf" >> "${filename}"
    fi
}

export -f experiment
export MAX_TIME
export ERRORS
export ALGORITHM
export HEURISTIC
export EXTRA

#Use GNU Parallel to run instances in parallel
parallel experiment ::: ${INSTANCE_LIST}

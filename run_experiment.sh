#!/bin/bash

MAX_TIME="30m"
ERRORS='errors.out'
INSTANCE_PATH='./instances'
INSTANCE_LIST=$(find ${INSTANCE_PATH} -type f -name '*.ideal')

ALGORITHM=$1
HEURISTIC=$2

function experiment {
    filename=$(basename -- $1)
    filename="${filename%.*}"
    filename="${filename}.test"
    timeout "${MAX_TIME}" sage experiment.sage $1 "${ALGORITHM}" "${HEURISTIC}" > "${filename}" 2>> "${ERRORS}"
    if [ $? -eq 124 ]
    then
        echo 'inf inf inf inf inf inf' >> "${filename}"
    fi
}

export -f experiment
export MAX_TIME
export ERRORS
export ALGORITHM
export HEURISTIC

#Use GNU Parallel to run instances in parallel
parallel experiment ::: ${INSTANCE_LIST}

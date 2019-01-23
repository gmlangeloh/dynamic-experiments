#!/bin/bash

#Usage example: run_summarize.sh cp-hilbert

EXPERIMENT_NAME=$1

#Put instance names in experiment files if missing
python fix_infs.py ${EXPERIMENT_NAME}

#Delete empty experiment result files
find ${EXPERIMENT_NAME} -type f -size 0 -delete

#Organize results in a single table
table="${EXPERIMENT_NAME}/${EXPERIMENT_NAME}.out"
cat ${EXPERIMENT_NAME}/*.test > ${table}

#Group results by number of variables
sage summarize.sage ${table} > results/"summary-${EXPERIMENT_NAME}.out"

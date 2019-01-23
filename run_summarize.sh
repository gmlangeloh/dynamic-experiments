#!/bin/bash

#Usage example: run_summarize.sh cp-hilbert

function generate_summary(){
  EXPERIMENT_NAME=$1
  echo -n "Generating summary for ${EXPERIMENT_NAME}..."
  #Put instance names in experiment files if missing
  python fix_infs.py "${EXPERIMENT_NAME}/"

  #Delete empty experiment result files
  find ${EXPERIMENT_NAME} -type f -size 0 -delete

  #Organize results in a single table
  table="${EXPERIMENT_NAME}/${EXPERIMENT_NAME}.out"
  cat ${EXPERIMENT_NAME}/*.test > "${table}"

  #Group results by number of variables
  sage summarize.sage ${table} > "results/summary-${EXPERIMENT_NAME}.out"
  echo "done"
  echo
}

ALGORITHMS="static
cp-hilbert
cp-betti
cp-mixed
caboara-hilbert
caboara-betti
caboara-mixed
gs-hilbert
simplex-hilbert
perturb-hilbert
perturb-random-hilbert
random-hilbert
random-betti
"

for algorithm in $ALGORITHMS
do
  generate_summary "$algorithm"
done

#Make LaTeX tables
Rscript generate_tables.r

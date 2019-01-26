#!/bin/bash

./run_experiment.sh random hilbert
mkdir random-hilbert
mv *.test random-hilbert

./run_experiment.sh random betti
mkdir random-betti
mv *.test random-betti

./run_experiment.sh random mixed
mkdir random-mixed
mv *.test random-mixed

./run_experiment.sh perturbation hilbert
mkdir perturb-hilbert
mv *.test perturb-hilbert

./run_experiment.sh perturbation betti
mkdir perturb-betti
mv *.test perturb-betti

./run_experiment.sh perturbation mixed
mkdir perturb-mixed
mv *.test perturb-mixed

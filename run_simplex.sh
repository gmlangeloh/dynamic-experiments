#!/bin/bash

./run_experiment.sh simplex hilbert
mkdir simplex-hilbert
mv *.test simplex-hilbert

./run_experiment.sh simplex betti
mkdir simplex-betti
mv *.test simplex-betti

./run_experiment.sh simplex mixed
mkdir simplex-mixed
mv *.test simplex-mixed

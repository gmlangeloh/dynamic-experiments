#!/bin/bash

./run_experiment.sh simplex-then-cp hilbert
mkdir simplex-then-cp-hilbert
mv *.test simplex-then-cp-hilbert

./run_experiment.sh simplex-then-cp betti
mkdir simplex-then-cp-betti
mv *.test simplex-then-cp-betti

./run_experiment.sh simplex-then-cp mixed
mkdir simplex-then-cp-mixed
mv *.test simplex-then-cp-mixed

./run_experiment.sh gs-then-cp hilbert
mkdir gs-then-cp-hilbert
mv *.test gs-then-cp-hilbert

./run_experiment.sh gs-then-cp betti
mkdir gs-then-cp-betti
mv *.test gs-then-cp-betti

./run_experiment.sh gs-then-cp mixed
mkdir gs-then-cp-mixed
mv *.test gs-then-cp-mixed

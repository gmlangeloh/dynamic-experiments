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

,/run_experiment.sh perturbation mixed
mkdir perturb-mixed
mv *.test perturb-mixed

./run_experiment.sh simplex hilbert
mkdir simplex-hilbert
mv *.test simplex-hilbert

./run_experiment.sh simplex betti
mkdir simplex-betti
mv *.test simplex-betti

,/run_experiment.sh simplex mixed
mkdir simplex-mixed
mv *.test simplex-mixed

./run_experiment.sh gritzmann-sturmfels hilbert
mkdir gs-hilbert
mv *.test gs-hilbert

./run_experiment.sh gritzmann-sturmfels betti
mkdir gs-betti
mv *.test gs-betti

./run_experiment.sh gritzmann-sturmfels mixed
mkdir gs-mixed
mv *.test gs-mixed

./run_experiment.sh regrets hilbert
mkdir regrets-hilbert
mv *.test regrets-hilbert

./run_experiment.sh regrets betti
mkdir regrets-betti
mv *.test regrets-betti

,/run_experiment.sh regrets mixed
mkdir regrets-mixed
mv *.test regrets-mixed

./run_experiment.sh gs-then-cp hilbert
mkdir gs-then-cp-hilbert
mv *.test gs-then-cp-hilbert

./run_experiment.sh gs-then-cp betti
mkdir gs-then-cp-betti
mv *.test gs-then-cp-betti

,/run_experiment.sh gs-then-cp mixed
mkdir gs-then-cp-mixed
mv *.test gs-then-cp-mixed

./run_experiment.sh simplex-then-cp hilbert
mkdir simplex-then-cp-hilbert
mv *.test simplex-then-cp-hilbert

./run_experiment.sh simplex-then-cp betti
mkdir simplex-then-cp-betti
mv *.test simplex-then-cp-betti

,/run_experiment.sh simplex-then-cp mixed
mkdir simplex-then-cp-mixed
mv *.test simplex-then-cp-mixed

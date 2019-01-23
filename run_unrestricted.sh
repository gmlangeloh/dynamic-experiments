#!/bin/bash

./run_experiment.sh gritzmann-sturmfels hilbert
mkdir gs-hilbert
mv *.test gs-hilbert

./run_experiment.sh gritzmann-sturmfels betti
mkdir gs-betti
mv *.test gs-betti

./run_experiment.sh gritzmann-sturmfels mixed
mkdir gs-mixed
mv *.test gs-mixed

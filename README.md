README
======

Repository for prototype implementations of dynamic Buchberger algorithms based
on the implementation of
[Caboara and Perry](http://www.math.usm.edu/perry/Research/dynamic_gb.pyx).
Includes additional *unrestricted* dynamic algorithms.

Everything requires SageMath (tested with SageMath 8.3).

Calling from SageMath
---------------------

Load the code with `load(buchberger.pyx)` and call the `dynamic_gb` function.

Running experiments
-------------------

Requires GNU Parallel. To run an experiment, run:

`./run_experiment.sh <algorithm> <heuristic>`

where <algorithm> is one of:
  - static
  - caboara-perry
  - random
  - caboara
  - perturbation
  - gritzmann-sturmfels
  - simplex
  - gs-then-cp
and <heuristic> is one of:
  - hilbert
  - betti
  - mixed

Making tables
-------------

Requires R/Rscript and Python. To make all tables, run:

`./run_summarize.sh`

Result files (output of `run_experiment.sh`) are expected to be in the following
directories:
  - static
  - cp-hilbert
  - cp-betti
  - cp-mixed
  - caboara-hilbert
  - caboara-betti
  - caboara-mixed
  - gs-hilbert
  - gs-betti
  - gs-mixed
  - random-hilbert
  - random-betti
  - random-mixed
  - perturb-hilbert
  - perturb-betti
  - perturb-mixed
  - simplex-hilbert
  - simplex-betti
  - simplex-mixed
  - regrets-hilbert
  - regrets-betti
  - regrets-mixed
  - gs-then-cp-hilbert
  - gs-then-cp-betti
  - gs-then-cp-mixed

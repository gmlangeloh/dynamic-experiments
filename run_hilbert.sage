'''
Tests the effects of making an initial analysis using the Hilbert heuristic.

If min_hilbert is passed as command line argument, outputs the basis sizes for
all orderings of minimal Hilbert degree instead.
'''

import sys

load("experiment.sage")
load("minkowski.sage")

if len(sys.argv) == 2 and sys.argv[1] == "min_hilbert":
    run_all_parallel('./instances/*.ideal', mindeg_hilbert_data)
else:
    run_all_parallel('./instances/*.ideal', general_hilbert_data)

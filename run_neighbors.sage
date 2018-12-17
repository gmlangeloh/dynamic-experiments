'''
Tests a quick neighbor count of some random initial point
of a Minkowski sum.

Precisely, we count how many changes of a single summand of the
initial vertex is an actual neighbor of it.
'''

load("experiment.sage")
load("feasible_neighbors.sage")

run_all_parallel('./instances/*.ideal', neighbor_data)

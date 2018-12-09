'''
Tests a simple version of the Betti heuristic as initial analysis.
'''

load("buchberger_graph.sage")
load("experiment.sage")

run_all_parallel('./instances/*.ideal', betti_heuristic_data)

'''
Computes Betti numbers with Bayer's criterion.
'''

load("buchberger_graph.sage")
load("experiment.sage")

run_all_parallel("./instances/*.ideal", bayer_data)

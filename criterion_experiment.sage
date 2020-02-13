'''
Runs some experiments on the divisibility criterion (Algorithm 1) of "A
new divisibility criterion to identify non-leading terms" by Mitchell
and Perry.

Applies the divisibility criterion to the local search algorithm.
'''

import sys

load("benchmarks.sage")
load("dynamicgb.pyx")

instances = [
    "cyclicn4",
    "cyclicnh4",
    "cyclicn5",
    "cyclicnh5",
    "cyclicn6",
    "cyclicnh6",
    "katsuran4",
    "katsuranh4",
    "katsuran5",
    "katsuranh5",
    "trinks",
    "econ4",
    "econh4",
    "econ5",
    "econh5"
]

extension = ".ideal"
directory = "./instances/"

criteria = [
    "newton",
    "perry1",
    "perry2"
]

for instance in instances:
    fullname = directory + instance + extension
    benchmark = Benchmark(fullname)
    for criterion in criteria:
        result = dynamic_gb(benchmark.ideal.gens(), algorithm="localsearch", return_stats = True, seed = 0, reducer = "classical", lscriterion = criterion, print_criterion=True)
        out = result[-1]
        print(instance, end=" ")
        print(criterion, end=" ")
        print(out)
        sys.stdout.flush()

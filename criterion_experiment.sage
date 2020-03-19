'''
Runs some experiments on the divisibility criterion (Algorithm 1) of "A
new divisibility criterion to identify non-leading terms" by Mitchell
and Perry.

Applies the divisibility criterion to the local search algorithm.
'''

import sys

load("benchmarks.sage")
load("new_reader.sage")
load("dynamicgb.pyx")

instances = [ [
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
], [
    "buch85",
    "butcher8",
    "eco8",
    "kotsireas",
    "noon3",
    "noon5",
    "s9_1"
] ]

extensions = [ ".ideal", ".txt" ]
directories = [ "./instances/", "./instances2/" ]

criteria = [
    "newton",
    "perry1",
    "perry2"
]

if len(sys.argv) > 1:
    tenure = int(sys.argv[1])
else:
    tenure = 0

for option in [1, 2]:
    extension = extensions[option]
    directory = directories[option]
    instance_list = instances[option]
    for instance in instance_list:
        fullname = directory + instance + extension
        if option == 1:
            benchmark = Benchmark(fullname)
            generators = benchmark.ideal.gens()
        elif option == 2:
            generators = read_ideal(fullname)
        for criterion in criteria:
            result = dynamic_gb(generators, algorithm="localsearch",
                                return_stats = True, seed = 0,
                                reducer = "classical", lscriterion = criterion,
                                print_criterion=True, taboo_tenure = tenure)
            out = result[-1]
            print(instance, end=" ")
            print(criterion, end=" ")
            print(out)
            sys.stdout.flush()

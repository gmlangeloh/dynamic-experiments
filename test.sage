'''
Test some dynamic GB algorithms.
'''

import sys

load("benchmarks.sage")
load("buchberger.pyx")

#Run the instance named sys.argv[1], if available. Else run a list of instances
if len(sys.argv) > 1:
  instances = [ sys.argv[1] ]
else:
  #These instances were chosen because they are small but non-trivial.
  #This means they are easy to test whenever changes are made to the codebase.
  instances = [ "cyclicn4",
                "cyclicnh4",
                "cyclicn5",
                "cyclicnh5",
                "cyclicn6",
                "cyclicnh6",
                "katsuran4",
                "katsuranh4",
                "katsuran5",
                "katsuranh5",
                "katsuran6",
                "katsuranh6" ]

if len(sys.argv) > 2:
  algorithms = [ sys.argv[2] ]
else:
  algorithms = [ "static",
                 #"caboara",
                 "caboara-perry",
                 "perturbation" ]

prefix = "./instances/"
suffix = ".ideal"

def instance_path(instance):
  return prefix + instance + suffix

#I am keeping here only the algorithms that seem promising

print("instance algorithm time dynamic heuristic queue reduction polynomials monomials degree sreductions zeroreductions")

for algorithm in algorithms:
  for instance in instances:
    benchmark = Benchmark(instance_path(instance))
    print(instance, end=" ")
    _ = dynamic_gb(benchmark.ideal.gens(), algorithm=algorithm, print_results=True)

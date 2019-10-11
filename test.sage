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
  basic_instances = [ "cyclicn4",
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
                      "katsuranh6"]
  #These instances are a bit larger, but should still complete in a few hours even
  #with the 30 repetitions for each algorithm.
  additional_instances = [ "cyclicn7",
                           "cyclicnh7",
                           "katsuran7",
                           "katsuranh7",
                           "butcher8",
                           "econ4",
                           "econh4",
                           "econ5",
                           "econh5",
                           "econ6",
                           "econh6",
                           "econ7",
                           "econh7",
                           "f633",
                           "f633h",
                           "virasoro",
                           "noon7"]
  instances = basic_instances + additional_instances

if len(sys.argv) > 2:
  algorithms = [ sys.argv[2] ]
else:
  #I am keeping here only the algorithms that seem promising
  algorithms = [ "static",
                 #"caboara",
                 "caboara-perry",
                 "perturbation",
                 "random" ]

prefix = "./instances/"
suffix = ".ideal"

def instance_path(instance):
  return prefix + instance + suffix

print("instance rep algorithm time dynamic heuristic queue reduction polynomials monomials degree sreductions zeroreductions")

for algorithm in algorithms:
  for instance in instances:
    benchmark = Benchmark(instance_path(instance))
    for repetition in range(30):
      print(instance, end=" ")
      print(repetition, end=" ")
      _ = dynamic_gb(benchmark.ideal.gens(), algorithm=algorithm, print_results=True)

'''
Test some dynamic GB algorithms.
'''

import io
import random
import sys
from contextlib import redirect_stdout
from multiprocessing import Pool, Lock

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
                 "caboara-perry",
                 "perturbation",
                 "random" ]

prefix = "./instances/"
suffix = ".ideal"

def instance_path(instance):
  return prefix + instance + suffix

def run_algorithm(algorithm, instance, repetition):

  #Run the experiment
  benchmark = Benchmark(instance_path(instance))
  f = io.StringIO()
  with redirect_stdout(f):
    _ = dynamic_gb(benchmark.ideal.gens(), algorithm=algorithm, print_results=True)
  out = f.getvalue()

  #Print correctly in stdout
  lock.acquire()
  try:
    print(instance, end=" ")
    print(repetition, end=" ")
    print(out, end="")
    sys.stdout.flush()
  finally:
    lock.release()

def init(l):
  global lock
  lock = l

print("instance rep algorithm time dynamic heuristic queue reduction polynomials monomials degree sreductions zeroreductions")

#Prepare (shuffled) list of experiments
lock = Lock()
experiments = []
for algorithm in algorithms:
  for instance in instances:
    for repetition in range(30):
      experiments.append((algorithm, instance, repetition))
random.shuffle(experiments)

with Pool(initializer=init, initargs=(lock,), processes=4) as pool:
  for experiment in experiments:
    pool.apply_async(run_algorithm, args = experiment)
  pool.close()
  pool.join()

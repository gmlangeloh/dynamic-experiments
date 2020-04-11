'''
Test some dynamic GB algorithms.
'''

import io
import random
import sys
from contextlib import redirect_stdout
from multiprocessing import Pool, Lock
from threading import Thread

import select_instances
load("benchmarks.sage")
load("dynamicgb.pyx")

#Run the instance named sys.argv[1], if available. Else run a list of instances
basic_instances_only = False
char0_only = False
easy = False
medium = False
full = False
single_instance = False
if len(sys.argv) > 1:
  if sys.argv[1] == "basic":
    basic_instances_only = True
  elif sys.argv[1] == "0":
    char0_only = True
  elif sys.argv[1] == "easy":
    easy = True
  elif sys.argv[1] == "medium":
    medium = True
  elif sys.argv[1] == "full":
    full = True
  else:
    single_instance = True
    instances = [ sys.argv[1] ]
if not single_instance:
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
                           "butcher8",
                           "fateman",
                           "fatemanh",
                           "noon7"]
  even_more_instances = [ "buch85.txt",
                          "butcher8.txt",
                          "eco8.txt",
                          "kotsireas.txt",
                          "noon3.txt",
                          "noon5.txt",
                          "s9_1.txt"]
  char0_instances = [
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
    "katsuran6",
    "katsuranh6",
    "cyclicn7",
    "cyclicnh7",
    "katsuran7",
    "katsuranh7",
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
    "butcher8",
    "fateman",
    "fatemanh",
    "noon7"
  ]
  if char0_only:
    instances = char0_instances
  elif basic_instances_only:
    instances = basic_instances
  else:
    instances = basic_instances + additional_instances

if easy or medium or full:
  easy, medium, toohard = select_instances.instance_selection()
  if easy:
    instances = easy
  elif medium:
    instances = medium
  elif full:
    instances = easy + medium

if len(sys.argv) > 2:
  algorithms = [ sys.argv[2] ]
else:
  #I am keeping here only the algorithms that seem promising
  algorithms = [ "static",
                 "caboara-perry",
                 "perturbation",
                 "random" ]
                 #"localsearch" ]

prefix = "./instances/"
suffix = ".ideal"

def instance_path(instance, coefs = False):
  if coefs:
    return "./instances-char0/" + instance + suffix
  return prefix + instance + suffix

def run_algorithm(experiment):

  algorithm, instance, reducer, repetition, coefs = experiment
  #Run the experiment
  timeout = 3600 #Use 1 hour timeout
  benchmark = Benchmark(instance_path(instance, coefs))
  result = dynamic_gb(benchmark.ideal.gens(), algorithm=algorithm, \
                      return_stats=True, seed=repetition, reducer=reducer, \
                      timeout=timeout, print_coefficients=coefs)
  #out = f.getvalue()
  out = result[-1]

  #Print correctly in stdout
  lock.acquire()
  try:
    print(instance, end=" ")
    print(reducer, end=" ")
    print(repetition, end=" ")
    print(out)
    sys.stdout.flush()
  finally:
    lock.release()

def init(l):
  global lock
  lock = l

def print_header(coefs = False):
  print("instance reducer rep algorithm time dynamic heuristic queue reduction polynomials monomials degree sreductions zeroreductions", end="")
  if coefs:
    print(" totalcoefs avgcoefs maxcoefs")
  else:
    print()

print_header()
#prepare list of experiments
lock = Lock()
experiments = []
reducers = [ 'classical'] if char0_only else [ 'classical', 'F4' ]

for instance in instances:
  for algorithm in algorithms:
    for reducer in reducers:
      if algorithm in ["random", "perturbation"]:
        for repetition in range(30):
          experiments.append((algorithm, instance, reducer, repetition, char0_only))
      else:
        experiments.append((algorithm, instance, reducer, 0, char0_only))

with Pool(initializer=init, initargs=(lock,), processes=4) as pool:
  pool.map(run_algorithm, experiments)

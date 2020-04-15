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
load("new_reader.sage")
load("dynamicgb.pyx")

#Run the instance named sys.argv[1], if available. Else run a list of instances
basic_instances_only = False
char0_only = False
easy = False
medium = False
full = False
extra = False
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
  elif sys.argv[1] == "extra":
    extra = True
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
  even_more_instances = [ "buch85",
                          "eco8",
                          "kotsireas",
                          "noon3",
                          "noon5",
                          "s9_1"]
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
  elif extra:
    instances = even_more_instances
  else:
    instances = basic_instances + additional_instances

if easy or medium or full:
  easy, medium, toohard = select_instances.instance_selection()
  if easy:
    instances = [ i for i in easy if i not in instances ]
  elif medium:
    instances = [ i for i in medium if i not in instances ]
  elif full:
    instances = [i for i in easy + medium if i not in instances ]

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

def instance_path(instance, coefs = False, extra = False):
  if coefs:
    return "./instances-char0/" + instance + suffix
  if extra:
    return "./instances2/" + instance + ".txt"
  return prefix + instance + suffix

def run_algorithm(experiment):

  algorithm, instance, reducer, repetition, coefs, extra = experiment
  #Run the experiment
  timeout = 3600 #Use 1 hour timeout
  if extra: #Running extra instances which are written in a different format
    generators = read_ideal(instance_path(instance, coefs, extra))
  else:
    benchmark = Benchmark(instance_path(instance, coefs, extra))
    generators = benchmark.ideal.gens()
  result = dynamic_gb(generators, algorithm=algorithm, \
                      return_stats=True, seed=repetition, reducer=reducer, \
                      timeout=timeout, print_coefficients=coefs)
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
          experiments.append((algorithm, instance, reducer, repetition, char0_only, extra))
      else:
        experiments.append((algorithm, instance, reducer, 0, char0_only, extra))
random.shuffle(experiments) #This is useful for better load balancing in the pool

with Pool(initializer=init, initargs=(lock,)) as pool:
  pool.map(run_algorithm, experiments)

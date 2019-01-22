'''
Generic tools to run experiments in Sage.
'''

import glob
import multiprocessing
import sys
from functools import partial
from multiprocessing.pool import Pool, ThreadPool

timeout = 30 * 60 #30 minutes timeout

def function_with_timeout(f, *args):
  p = ThreadPool(1)
  res = p.apply_async(f, args)
  try:
    out = res.get(timeout)
    return out
  except multiprocessing.TimeoutError:
    p.terminate()
    #Print same amount of columns as in the experiments
    print "inf inf inf inf inf inf"
  except Exception as e:
    sys.stderr.write(str(e) + "\n")
  return None

def run_all_parallel(glob_pattern, experiment_function):
  instances = glob.glob(glob_pattern)
  pool = Pool()
  timed_experiment = partial(function_with_timeout, experiment_function)
  pool.map(timed_experiment, instances)

def valid_instance(benchmark):
  return benchmark.ideal.ring().ngens() <= 8

def instance_name(instance_path):
  return instance_path.split("/")[-1].split(".")[0]

def run_static(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    #f = open(name+".out", "w")
    #sys.stdout = f
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", static=True, \
                       print_results=True)
    sys.stdout.flush()

def run_caboara_perry(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    #f = open(name+".out", "w")
    #sys.stdout = f
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", print_results=True, \
                       heuristic=global_heuristic)
    sys.stdout.flush()

def run_caboara(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    #f = open(name+".out", "w")
    #sys.stdout = f
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", print_results=True, \
                       heuristic=global_heuristic, use_boundary_vectors=False,\
                       use_disjoint_cones=False)
    sys.stdout.flush()

def run_random(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    #f = open(name+".out", "w")
    #sys.stdout = f
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", print_results=True, \
                       heuristic=global_heuristic, random=True)
    sys.stdout.flush()

def run_perturbation(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", print_results=True, \
                       heuristic=global_heuristic, perturbation=True)
    sys.stdout.flush()

def run_gritzmann_sturmfels(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", print_results=True, \
                       heuristic=global_heuristic, unrestricted=True)
    sys.stdout.flush()

instance_glob = './instances/*.ideal'
global_heuristic = 'hilbert'

load("benchmarks.sage")
if not valid_instance(Benchmark(sys.argv[1])):
  sys.exit(0)

sys.stderr.write("Starting: " + sys.argv[1] + "\n")
sys.stderr.flush()
if len(sys.argv) > 3:
  if sys.argv[3] in ['hilbert', 'betti', 'mixed']:
    global_heuristic = sys.argv[3]

if len(sys.argv) > 2:
  load("buchberger.pyx")
  if sys.argv[2] == 'static':
    run_static(sys.argv[1])
  elif sys.argv[2] == 'caboara-perry':
    run_caboara_perry(sys.argv[1])
  elif sys.argv[2] == 'random':
    run_random(sys.argv[1])
  elif sys.argv[2] == 'caboara':
    run_caboara(sys.argv[1])
  elif sys.argv[2] == 'perturbation':
    run_perturbation(sys.argv[1])
  elif sys.argv[2] == 'gritzmann-sturmfels':
    run_gritzmann_sturmfels(sys.argv[1])

sys.stderr.write("Finished" + sys.argv[1] + "\n")
sys.stderr.flush()

'''
Generic tools to run experiments in Sage.
'''

import glob
import multiprocessing
import sys
from functools import partial
from multiprocessing.pool import Pool

timeout = 30 * 60 #30 minutes timeout

def function_with_timeout(f, args):
  p = Pool(1)
  res = p.apply_async(f, args)
  try:
    out = res.get(timeout)
    return out
  except multiprocessing.TimeoutError:
    p.terminate()
    #Print same amount of columns as in the experiments
    print "inf inf inf inf inf inf"
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
    print instance_name(instance_path),
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", static=True, \
                       print_results=True)

def run_caboara_perry(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    print instance_name(instance_path),
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", print_results=True \
                       heuristic=heuristic)

instance_glob = './instances/*.ideal'
global_heuristic = 'hilbert'

if len(sys.argv) > 2:
  if sys.argv[2] in ['hilbert', 'betti', 'mixed']
    global_heuristic = sys.argv[2]

if len(sys.argv) > 1:
  load("buchberger.pyx")
  if sys.argv[1] == 'static':
    run_all_parallel(instance_glob, run_static)
  elif sys.argv[1] == 'caboara-perry':
    run_all_parallel(instance_glob, run_caboara_perry)

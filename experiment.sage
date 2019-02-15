'''
Generic tools to run experiments in Sage.
'''

import glob
import sys

def valid_instance(benchmark):
  return benchmark.ideal.ring().ngens() <= 8

def instance_name(instance_path):
  return instance_path.split("/")[-1].split(".")[0]

def run_static(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", static=True, \
                       print_results=True)
    sys.stdout.flush()

def run_caboara_perry(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", print_results=True, \
                       heuristic=global_heuristic)
    sys.stdout.flush()

def run_caboara(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", print_results=True, \
                       heuristic=global_heuristic, use_boundary_vectors=False,\
                       use_disjoint_cones=False)
    sys.stdout.flush()

def run_random(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
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
                       heuristic=global_heuristic, perturbation=True, \
                       initial_ordering=init_order)
    sys.stdout.flush()

def run_gritzmann_sturmfels(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", print_results=True, \
                       heuristic=global_heuristic, unrestricted=True)
    sys.stdout.flush()

def run_simplex(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", print_results=True, \
                       heuristic=global_heuristic, simplex=True, \
                       initial_ordering=init_order)
    sys.stdout.flush()

def run_regrets(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy="sugar", print_results=True, \
                       heuristic=global_heuristic, reinsert=True, \
                       use_disjoint_cones=False, use_boundary_vectors=False)
    sys.stdout.flush()

def run_gs_then_cp(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy='sugar', print_results=True, \
                       heuristic=global_heuristic, unrestricted=True, \
                       max_calls=len(B.ideal.gens()))
    sys.stdout.flush()

def run_simplex_then_cp(instance_path):
  B = Benchmark(instance_path)
  if valid_instance(B):
    name = instance_name(instance_path)
    print name,
    dummy = dynamic_gb(B.ideal.gens(), strategy='sugar', print_results=True, \
                       heuristic=global_heuristic, simplex=True, \
                       max_calls=20)
    sys.stdout.flush()

instance_glob = './instances/*.ideal'
global_heuristic = 'hilbert'
init_order = 'grevlex'

load("benchmarks.sage")
if not valid_instance(Benchmark(sys.argv[1])):
  sys.exit(0)

sys.stderr.write("Starting: " + sys.argv[1] + "\n")
sys.stderr.flush()

if len(sys.argv) > 3:
  if sys.argv[3] in ['hilbert', 'betti', 'mixed']:
    global_heuristic = sys.argv[3]

if len(sys.argv) > 4:
  if sys.argv[4] in ['grevlex', 'init_random']:
    init_order = sys.argv[4]

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
  elif sys.argv[2] == 'simplex':
    run_simplex(sys.argv[1])
  elif sys.argv[2] == 'regrets':
    run_regrets(sys.argv[1])
  elif sys.argv[2] == 'gs-then-cp':
    run_gs_then_cp(sys.argv[1])
  elif sys.argv[2] == 'simplex-then-cp':
    run_simplex_then_cp(sys.argv[1])

sys.stderr.write("Finished" + sys.argv[1] + "\n")
sys.stderr.flush()

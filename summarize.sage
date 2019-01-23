import glob
import sys

if len(sys.argv) <= 1:
  print "Results file missing."
  exit(1)

load("benchmarks.sage")

path_to_instances = "./instances/*.ideal"
instances = glob.glob(path_to_instances)

#Build a dictionary { instances -> nvars }
num_vars = {}
for instance in instances:
  B = Benchmark(instance)
  instance_name = instance.split("/")[-1].split(".")[0]
  num_vars[instance_name] = B.ideal.ring().ngens()

#Summarize results per instance size:
class NVarResults:

  def __init__(self, n):
    self.n = n
    self.time = 0.0
    self.overhead = 0.0
    self.polynomials = 0
    self.monomials = 0
    self.max_degree = 0
    self.reductions = 0

    self.instances = []
    self.timeouts = 0

  def update(self, instance_results):
    results = instance_results.split()
    instance_name = results[0]
    if results[1] != 'inf':
      self.time += float(results[1])
      self.overhead += float(results[2])
      self.polynomials += int(results[3])
      self.monomials += int(results[4])
      self.max_degree += int(results[5])
      self.reductions += int(results[6])
    else:
      self.timeouts += 1

    self.instances.append(instance_name)

  def _average(self):
    i = float(len(self.instances))
    return self.time / i, self.overhead / i, self.polynomials / i, \
        self.monomials / i, self.max_degree / i, self.reductions / i

  def report(self):
    time, overhead, polys, monomials, deg, reductions = self._average()
    print self.n, len(self.instances), self.timeouts, time, overhead, polys, monomials, deg, reductions

def header():
  print "Variables Instances Timeouts time overhead #G monomials maxdeg reductions"

min_nvars = 2
max_nvars = 8

results = {}
for n in xrange(min_nvars, max_nvars + 1):
  results[n] = NVarResults(n)

with open(sys.argv[1], "r") as f:
  for line in f:
    instance_name = line.split()[0]
    n = num_vars[instance_name]
    results[n].update(line)

header()
for n in xrange(min_nvars, max_nvars + 1):
  results[n].report()

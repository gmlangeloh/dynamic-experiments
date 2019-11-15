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

    self.time_timeouts = 0.0

    self.instances = []
    self.timeouts = 0

  def update(self, instance_results):
    results = instance_results.split()
    instance_name = results[0]
    if len(results) < 2:
      sys.stderr.write("problems with " + instance_name + "\n")
    if results[1] != 'inf':
      self.time += float(results[1])
      self.overhead += float(results[2])
      self.polynomials += int(results[3])
      self.monomials += int(results[4])
      self.max_degree += int(results[5])
      self.reductions += int(results[6])

      self.time_timeouts += float(results[1])
    else:
      self.timeouts += 1
      self.time_timeouts += 30.0 * 60.0 #The timeout time

    self.instances.append(instance_name)

  def _average(self):
    i = float(len(self.instances))
    j = float(len(self.instances) - self.timeouts)
    if len(self.instances) == self.timeouts:
      Inf = float('inf')
      sys.stderr.write("No instances of size: " + str(self.n) + "\n")
      return self.time_timeouts / i, Inf, Inf, Inf, Inf, Inf, Inf
    return self.time_timeouts / i , self.time / j, self.overhead / j, \
      self.polynomials / j, self.monomials / j, self.max_degree / j, \
      self.reductions / j

  def report(self):
    t1, t2, overhead, polys, monomials, deg, reductions = self._average()
    print self.n, self.timeouts, t1, t2, overhead, polys, monomials, deg, reductions

def header():
  print "Variables Timeouts t1 t2 overhead polys monomials maxdeg reductions"

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


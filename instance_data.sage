import glob
import sys

load("benchmarks.sage")

path_to_instances = "./instances/*.ideal"
instances = glob.glob(path_to_instances)
instances.sort()

def instance_name(path):
  return path.split("/")[-1].split(".")[0]

def instance_header():
  print 'name n m k T deg'

def summary_header():
  print 'n m T deg'

def valid_benchmark(B):
  return B.ideal.ring().ngens() <= 8

class NVarData:

  def __init__(self, n):
    self.n = n
    self.m = 0
    self.mon = 0
    self.deg = 0

    self.num_instances = 0

  def _average(self):
    i = float(self.num_instances)
    return self.m / i, self.mon / i, self.deg / i

  def update(self, B):
    n, m, k, mon, deg = B.data()
    assert n == self.n, "Wrong number of variables in this group"
    self.m += m
    self.mon += mon
    self.deg += deg

    self.num_instances += 1

  def report(self):
    polys, mons, degs = self._average()
    print self.n, polys, mons, degs

min_nvars = 2
max_nvars = 8

nvars = {}
for n in xrange(min_nvars, max_nvars + 1):
  nvars[n] = NVarData(n)

with open("instances.out", "w") as sys.stdout:
  #Print data about each instance
  for instance in instances:
    B = Benchmark(instance)
    if valid_benchmark(B):
      name = instance_name(instance)
      n, m, k, mon, deg = B.data()
      nvars[n].update(B)
      print name, n, m, k, mon, deg

with open("summary-instances.out", "w") as sys.stdout:
  #Print data about each instance group (averages, etc)
  for n in xrange(min_nvars, max_nvars + 1):
    nvars[n].report()

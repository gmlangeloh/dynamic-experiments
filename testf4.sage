'''
Run some F4 experiments.
'''

load("benchmarks.sage")
load("buchberger.pyx")

if len(sys.argv) > 1:
  instances = [ sys.argv[1] ]
else:
  instances = [ "cyclicn4",
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
                "katsuranh6" ]

prefix = "./instances/"
suffix = ".ideal"

def instance_path(instance):
  return prefix + instance + suffix

for instance in instances:
  benchmark = Benchmark(instance_path(instance))
  print(instance, end=" ")
  for reducer in ['classical', 'F4']:
    _ = dynamic_gb(benchmark.ideal.gens(), print_results=True, reducer=reducer)

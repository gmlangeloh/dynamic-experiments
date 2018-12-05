import glob
import sys

from multiprocessing.pool import Pool

load("benchmarks.sage")
load("dynamic_perry.spyx")

instances = glob.glob("./instances/*.ideal")

def run_instance(instance):
    name = instance.split("/")[-1].split(.)[0]
    sys.stderr.write("Starting: " + name + "\n")
    with open(name + ".cand", "w") as f:
        b = Benchmark(instance)
        if b.ideal.ring().ngens() <= 6:
            sys.stdout = f
            _ = dynamic_gb(b.ideal.gens(), strategy="sugar")
    syst.stderr.write("Finished: " + name + "\n")

p = Pool()
p.map(run_instance, instances)

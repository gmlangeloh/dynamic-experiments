from time import time
import glob
import sys
import multiprocessing

load("benchmarks.sage")

filenames = glob.glob("./instances/*.ideal")
exceptions = [ "./instances/cyclicnh9.ideal" ]

def maxdeg(G):
    return max([g.degree() for g in G])

def totalpolys(G):
    return len(G)

def totalmonoms(G):
    return sum([len(g.monomials()) for g in G])

def run_singular(filename):
    if filename in exceptions:
        return
    I = Benchmark(filename)
    t = time()
    G = I.ideal.groebner_basis()
    t = time() - t
    print(filename, t, totalpolys(G), totalmonoms(G), maxdeg(G))
    sys.stdout.flush()

for filename in filenames:
    p = multiprocessing.Process(target=run_singular, args=(filename,))
    p.start()
    p.join(3600)
    if p.is_alive():
        p.terminate()
        p.join()
        print(filename, "timed out at 3600s")

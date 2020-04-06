from time import time
import glob
import sys

load("benchmarks.sage")

filenames = glob.glob("./instances/*.ideal")
exceptions = []

def maxdeg(G):
    return max([g.degree() for g in G])

def totalpolys(G):
    return len(G)

def totalmonoms(G):
    return sum([len(g.monomials()) for g in G])

for f in filenames:
    if f not in exceptions:
        I = Benchmark(f)
        t = time()
        G = I.ideal.groebner_basis()
        t = time() - t
        print(f, t, totalpolys(G), totalmonoms(G), maxdeg(G))
        sys.stdout.flush()

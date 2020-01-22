import sys
load("benchmarks.sage")

def katsura_order(n):
    w = [2**(n-1)] * (n+1)
    w[n-1] = 1
    w[n] = 1
    return w

def basis_size_degree(G):
    size = len(G)
    degree = max([g.degree() for g in G])
    return size, degree

def katsura(n):
    path = "./instances/"
    name = "katsuran" + str(n)
    extension = ".ideal"
    filename = path + name + extension
    kat = Benchmark(filename)
    return kat

if len(sys.argv) >= 2:
    max_n = int(sys.argv[1])
else:
    max_n = 8

for n in range(3, max_n):
    kat = katsura(n)
    w = katsura_order(n)
    kat.change_order(w)
    G = kat.ideal.groebner_basis()
    size, degree = basis_size_degree(G)
    print(n, size, degree)

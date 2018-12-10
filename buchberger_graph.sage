'''
Experimentation with the Buchberger graph (see Miller and Sturmfels' 'Combinatorial Commutative Algebra' for a definition).
'''

import sys

load("benchmarks.sage")
load("minkowski.sage")

def is_edge(i, j, LMs):
    R = LMs[0].parent()
    l = R.monomial_lcm(LMs[i], LMs[j])
    for k in xrange(len(LMs)):
        if k != i and k != j and R.monomial_divides(LMs[k], l):
            if all( l.degree(v) == 0 or l.degree(v) > LMs[k].degree(v) \
                    for v in R.gens()):
                return False
    return True

def graph_edges(LMs):
    num_edges = 0
    for j in xrange(len(LMs)):
        for i in xrange(j):
            if is_edge(i, j, LMs):
                num_edges += 1
    return num_edges

def betti_heuristic_data(instance):
    name = instance.split("/")[-1].split(".")[0]
    I = Benchmark(instance).ideal
    if I.ring().ngens() > 8:
        return
    sys.stderr.write("Starting: " + name + "\n")
    Mink = minkowski(I)
    with open(name + ".betti", "w") as f:
        sys.stdout = f
        for v in Mink.vertex_generator():
            w = normal(v, Mink)
            R = PolynomialRing(I.ring().base_ring(), I.ring().gens(), order=create_order(w))
            G = [ R(g) for g in I.gens() ]
            LMs = [ g.lm() for g in G ]
            try:
                print graph_edges(LMs), len(R.ideal(G).groebner_basis())
            except:
                print graph_edges(LMs), 10000
    sys.stderr.write("Finished: " + name + "\n")

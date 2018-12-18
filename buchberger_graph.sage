'''
Experimentation with the Buchberger graph (see Miller and Sturmfels'
'Combinatorial Commutative Algebra' for a definition).

Also implements Bayer's graph criterion for minimal generating sets of Syzygy
modules.
'''

import sys

load("benchmarks.sage")
load("minkowski.sage")

def is_edge(i, j, LMs):
    R = LMs[0].parent()
    l = R.monomial_lcm(LMs[i], LMs[j])
    for k in xrange(len(LMs)):
        #TODO checking k != i and k != j is probably unnecessary
        #if k != i and k != j and R.monomial_divides(LMs[k], l):
        if R.monomial_divides(LMs[k], l):
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
            R = PolynomialRing(I.ring().base_ring(), I.ring().gens(), \
                               order=create_order(w))
            G = [ R(g) for g in I.gens() ]
            LMs = [ g.lm() for g in G ]
            try:
                print graph_edges(LMs), len(R.ideal(G).groebner_basis())
            except:
                print graph_edges(LMs), 10000
    sys.stderr.write("Finished: " + name + "\n")

def pair_graph(p, LMs, previously_eliminated):

    R = LMs[0].parent()
    k = len(LMs)
    a = p[0]
    b = p[1]
    l = R.monomial_lcm(LMs[a], LMs[b])
    vertices = [ i for i in xrange(k) if R.monomial_divides(LMs[i], l) ]
    edges = previously_eliminated
    for i in vertices:
        for j in vertices:
            if i != j and R.monomial_lcm(LMs[i], LMs[j]) != l:
                edges.append((i, j))

    return Graph(edges)

def bayer(LMs):
    '''
    Computes a minimal generating set of Syz(ideal(LMs))
    '''
    k = len(LMs)
    pairs = [ (i, j) for i in xrange(k) for j in xrange(i) ]
    eliminated = []
    minimal_gens = []
    for p in pairs:
        G = pair_graph(p, LMs, eliminated)
        try:
            if G.distance(p[0], p[1]) != Infinity:
                eliminated.append(p)
            else:
                minimal_gens.append(p)
        except ValueError:
            #p[0] or p[1] not in graph -> no edges incident to one of them
            #-> no path between p[0] and p[1]
            minimal_gens.append(p)
    return minimal_gens

def bayer_data(instance):
    name = instance.split("/")[-1].split(".")[0]
    I = Benchmark(instance).ideal
    if I.ring().ngens() > 7:
        return
    sys.stderr.write("Starting: " + name + "\n")
    Mink = minkowski(I)
    with open(name + ".bayer", "w") as f:
        sys.stdout = f
        for v in Mink.vertex_generator():
            w = normal(v, Mink)
            R = PolynomialRing(I.ring().base_ring(), I.ring().gens(), \
                               order=create_order(w))
            G = [ R(g) for g in I.gens() ]
            LMs = [ g.lm() for g in G ]
            betti = len(bayer(LMs))
            try:
                print I.ring().ngens(), betti, len(R.ideal(G).groebner_basis())
            except:
                print I.ring().ngens(), betti, 10000
    sys.stderr.write("Finished: " + name + "\n")

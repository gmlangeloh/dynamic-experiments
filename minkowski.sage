'''
Useful module to test the effects of making an initial analysis

We can print:
- best/worst Hilbert Functions, average
- compare with the ordering Perry's algorithm would choose in the same number of steps
- compare with grevlex as well
- compare size of final basis taking the chosen order statically
'''

import glob
import os
import sys

from multiprocessing.pool import Pool

load("benchmarks.sage")
load("dynamic_perry.spyx")

ZZ = IntegerRing()

def minkowski(I):
    n = I.ring().ngens()
    result = Polyhedron(rays=(-identity_matrix(n)).rows())
    for g in I.gens():
        result += g.newton_polytope()
    return result

def normal(v, P):
    inequalities = P.inequalities()
    indices = [ ieq.index() for ieq in v.incident() ]
    rays = [ -inequalities[i].A() for i in indices ]
    return list(sum(Cone(rays)))

def create_order(w):
    r"""
    Create term ordering acceptable to Singular, using integer weight vector ``w``.
    """
    # first we convert floats to integer
    # this is fine since we are using integer programming to find the ordering
    wZZ = [ZZ(each) for each in w]
    M = list()
    M.append(wZZ)

    # now fill in the lower rows of the matrix, to avoid potential ambiguity
    for i in xrange(len(w)-1):
        M.append([1 for k in xrange(i+1,len(w))] + [0 for k in xrange(i+1)])

    return TermOrder(matrix(M))

def evaluate(v, P, I):
    w = normal(v, P)
    R = PolynomialRing(I.ring().base_ring(), I.ring().gens(), order=create_order(w))
    G = [ R(g) for g in I.gens() ]
    LMs = [ g.lm() for g in G ]
    HP = R.ideal(LMs).hilbert_polynomial()
    return G, w, HP.degree(), HP.lc()

def evaluate_all(I):
    '''
    Prints information about all distinguishable orders for input ideal I
    '''
    Mink = minkowski(I)
    mindeg = float("inf")
    maxdeg = float("-inf")
    minw = []
    for v in Mink.vertex_generator():
        G, w, deg, coef = evaluate(v, Mink, I)
        LMs = [ g.lm() for g in G]
        #print LMs, w, deg, coef
        if deg < mindeg:
            mindeg = deg
            minw = w
        if deg > maxdeg:
            maxdeg = deg
    R = PolynomialRing(I.ring().base_ring(), I.ring().gens(), order=create_order(minw))
    G = [ R(g) for g in I.gens() ]
    print mindeg, maxdeg, len(R.ideal(G).groebner_basis()),

def evaluate_all_min(I):
    Mink = minkowski(I)
    mindeg = float("inf")
    minws = []
    for v in Mink.vertex_generator():
        G, w, deg, coef = evaluate(v, Mink, I)
        LMs = [ g.lm() for g in G ]
        if deg < mindeg:
            mindeg = deg
            minws = [w, coef]
        elif deg == mindeg:
            minws.append(w)
    for w, coef in minws:
        R = PolynomialRing(I.ring().base_ring(), I.ring().gens(), order=create_order(w))
        G = [ R(g) for g in I.gens() ]
        print coef, len(R.ideal(G).groebner_basis())

def evaluate_grevlex(I):
    '''
    Prints information about the grevlex order for I
    '''
    n = I.ring().ngens()
    w = [1] * n
    R = PolynomialRing(I.ring().base_ring(), I.ring().gens(), order=create_order(w))
    G = [ R(g) for g in I.gens() ]
    LMs = [ g.lm() for g in G ]
    HP = R.ideal(LMs).hilbert_polynomial()
    print HP.degree(), HP.lc(), len(R.ideal(G).groebner_basis()),

def evaluate_perry(I):
    '''
    Runs initial iterations of Perry's algorithm to obtain its `initial` ordering
    '''
    w = dynamic_gb(I.gens(), strategy="sugar", itmax=len(I.gens()))[0]
    R = PolynomialRing(I.ring().base_ring(), I.ring().gens(), order=create_order(w))
    G = [ R(g) for g in I.gens() ]
    LMs = [ g.lm() for g in G ]
    HP = R.ideal(LMs).hilbert_polynomial()
    print HP.degree(), HP.lc(), len(R.ideal(G).groebner_basis()),

def run_instance(instance):
    name = instance.split("/")[2].split(".")[0]
    sys.stderr.write("starting: " + name)
    b = Benchmark(instance)
    if b.ideal.ring().ngens() <= 8: #Try only relatively small instances
        print name,
        evaluate_grevlex(b.ideal)
        evaluate_perry(b.ideal)
        evaluate_all(b.ideal)
        print
    sys.stderr.write("finished: " + name)
    sys.stdout.flush()

def run_instance_min(instance):
    '''
    Run all orders with min Hilbert degree for instance
    '''
    name = instance.split("/")[2].split(".")[0]
    b = Benchmark(instance)
    if b.ideal.ring().ngens() <= 7: #Try only relatively small instances
        sys.stderr.write("starting: " + name)
        with open(name + ".min") as f:
            sys.stdout = f
            evaluate_all_min(b.ideal)
        sys.stderr.write("finished: " + name)

def run_all():
    instances = glob.glob('./instances/*.ideal')
    pool = Pool(processes=8)
    pool.map(run_instance, instances)

def run_all_min():
    instances = glob.glob('./instances/*.ideal')
    pool = Pool()
    pool.map(run_instance_min, instances)

run_all_min()

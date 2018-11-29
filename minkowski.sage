'''
Useful module to test the effects of making an initial analysis

We can print:
- best/worst Hilbert Functions, average
- compare with the ordering Perry's algorithm would choose in the same number of steps
- compare with grevlex as well
- compare size of final basis taking the chosen order statically
'''

from sage.rings.integer_ring import IntegerRing

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix

from sage.geometry.polyhedron.constructor import Polyhedron
from sage.geometry.cone import Cone

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder

from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular

import glob

ZZ = IntegerRing()

cpdef minkowski(I):
    cdef int n = I.ring().ngens()
    result = Polyhedron(rays=(-identity_matrix(n)).rows())
    for g in I.gens():
        result += g.newton_polytope()
    return result

cpdef normal(v, P):
    cdef tuple inequalities = P.inequalities()
    cdef list indices = [ ieq.index() for ieq in v.incident() ]
    cdef list rays = [ -inequalities[i].A() for i in indices ]
    return list(sum(Cone(rays)))

cpdef create_order(list w):
    r"""
    Create term ordering acceptable to Singular, using integer weight vector ``w``.
    """
    # first we convert floats to integer
    # this is fine since we are using integer programming to find the ordering
    cdef list wZZ = [ZZ(each) for each in w]
    cdef list M = list()
    M.append(wZZ)

    # now fill in the lower rows of the matrix, to avoid potential ambiguity
    cdef int i
    for i in xrange(len(w)-1):
        M.append([1 for k in xrange(i+1,len(w))] + [0 for k in xrange(i+1)])

    return TermOrder(matrix(M))

cpdef tuple evaluate(v, P, I):
    cdef list w = normal(v, P)
    cdef MPolynomialRing_libsingular R = PolynomialRing(I.ring().base_ring(), I.ring().gens(), order=create_order(w))
    cdef list G = [ R(g) for g in I.gens() ]
    cdef list LMs = [ g.lm() for g in G ]
    HP = R.ideal(LMs).hilbert_polynomial()
    return G, w, HP.degree(), HP.lc()

cpdef void evaluate_all(I):
    Mink = minkowski(I)
    cdef float mindeg = float("inf")
    cdef float maxdeg = float("-inf")
    cdef list minw = []
    for v in Mink.vertex_generator():
        G, w, deg, coef = evaluate(v, Mink, I)
        LMs = [ g.lm() for g in G]
        #print LMs, w, deg, coef
        if deg < mindeg:
            mindeg = deg
            minw = w
        if deg > maxdeg:
            maxdeg = deg
    cdef MPolynomialRing_libsingular R = PolynomialRing(I.ring().base_ring(), I.ring().gens(), order=create_order(minw))
    G = [ R(g) for g in I.gens() ]
    print mindeg, maxdeg, len(R.ideal(G).groebner_basis()),

cpdef void evaluate_grevlex(I):
    cdef int n = I.ring().ngens()
    cdef list w = [1] * n
    cdef MPolynomialRing_libsingular R = PolynomialRing(I.ring().base_ring(), I.ring().gens(), order=create_order(w))
    cdef list G = [ R(g) for g in I.gens() ]
    cdef list LMs = [ g.lm() for g in G ]
    HP = R.ideal(LMs).hilbert_polynomial()
    print HP.degree(), HP.lc(), len(R.ideal(G).groebner_basis()),

cpdef void evaluate_perry(I):
    pass

def run_all():
    instances = glob.glob('./instances/*.ideal')
    for instance in instances:
        name = instance.split("/")[2].split(".")[0]
        b = Benchmark(instance)
        if b.ideal.ring().ngens() <= 8: #Try only relatively small instances
            print name,
            evaluate_grevlex(b.ideal)
            evaluate_perry(b.ideal)
            evaluate_all(b.ideal)

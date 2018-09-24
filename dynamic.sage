'''
Simple implementation of Dynamic (unrestricted) Buchberger.
'''

from random import randint

def negative_orthant(n):
    '''
    The negative orthant of R^n as a polyhedron
    '''
    return Polyhedron(rays=(-identity_matrix(n)).rows())

class BasisElement:
    '''
    A BasisElement is a polynomial in a (partial) GB along with its affine
    Newton Polyhedron.
    '''
    def __init__(self, f):
        self._newton_polyhedron = None
        self._polynomial = f

    def polynomial(self):
        return self._polynomial

    def newton_polyhedron(self):
        '''
        The affine Newton polyhedron of a polynomial is the Minkowski sum of
        its Newton polytope with the negative orthant of R^n.

        Its vertices are in bijection with the possible leading monomials of
        the input polynomial.
        '''
        if self._newton_polyhedron is None:
            n = self.polynomial().parent().ngens()
            self._newton_polyhedron = self.polynomial().newton_polytope() + \
                negative_orthant(n)
        return self._newton_polyhedron

    def _normal_cone(self, v):
        rays = []
        for u in self.newton_polyhedron():
            if u != v:
                rays.append(vector(v) - vector(u))
        return Cone(rays)

    def normal_fan(self):
        '''
        The normal fan of this polynomial's Newton polyhedron as a list of
        normal cones of its vertices.
        '''
        fan = {}
        for v in self.newton_polyhedron().vertices():
            fan[tuple(v)] = self._normal_cone(v)
        return fan

    def change_order(self, w):
        new_order = TermOrder("wdegrevlex", w)
        new_ring = self._polynomial.parent().change_ring(order=new_order)
        self._polynomial = new_ring(self._polynomial)

class MonomialIdeal:
    '''
    A MonomialIdeal is used to keep an ideal with its associated Hilbert
    Function. MonomialIdeals can be compared with respect to their Hilbert
    Polynomials.
    '''
    def __init__(self, I, w):
        self._ideal = I
        self._weights = w
        self._hilbert = None

    def hilbert(self):
        if self._hilbert is None:
            self._hilbert = self._ideal.hilbert_polynomial()
        return self._hilbert

    def weights(self):
        return self._weights

    def __cmp__(self, other):
        '''
        self < other iff the degree of its Hilbert Polynomial is smaller or
        degrees are the same and leading coefficient is smaller.
        '''
        HP1 = self.hilbert()
        HP2 = other.hilbert()
        if HP1.degree() < HP2.degree():
            return -1
        elif HP1.degree() == HP2.degree():
            if HP1.lc() < HP2.lc():
                return -1
            elif HP1.lc() == HP2.lc():
                return 0
        return 1

class DynamicEngine:

    def __init__(self, R):
        n = R.ngens()
        self._order = TermOrder("wdegrevlex", [1] * n)
        self._n = n
        self._ring = R
        self._best_ideal = None
        self._call = 0

    def order(self):
        return _order

    def _random_vector(self):
        coords = [ randint(1, 100) for i in range(self._n)]
        return vector(coords)

    def _random_minkowski_vertex(self, G):
        w = self._random_vector()
        minkowski_vertex = vector([0] * self._n)
        vertex_decomposition = []
        for g in G:
            NP = g.newton_polyhedron()
            lp, var = NP.to_linear_program(return_variable = True)
            lp.set_objective(sum([ w[i] * var[i] for i in range(self._n)]))
            lp.solve()
            values = lp.get_values(var)
            vertex = [ values[i] for i in range(self._n) ]
            minkowski_vertex += vector(vertex)
            vertex_decomposition.append(vector(vertex))
        return list(w), minkowski_vertex, vertex_decomposition

    def _monomial_from_vertex(self, v):
        return prod([ self._ring.gens()[i]^v[i] for i in range(self._n) ])

    def _ideal_from_decomposition(self, decomposition):
        G = [ self._monomial_from_vertex(v) for v in decomposition ]
        return ideal(G)

    def change_order(self, w):
        self._order = TermOrder("wdegrevlex", w)

    #Unless I'm mistaken, we can also implement the restricted algorithm
    #easily working with Minkowski sums: it is enough to sum the previously
    #accepted vertex with the polytope of the new basis element!
    def next(self, G, iterations, period):
        self._call += 1
        if self._call % period != 1:
            return
        for i in range(iterations):
            w, v, v_decomposed = self._random_minkowski_vertex(G)
            I = MonomialIdeal(self._ideal_from_decomposition(v_decomposed), w)
            if self._best_ideal is None or I < self._best_ideal:
                self._best_ideal = I
        best_order = self._best_ideal.weights()
        self.change_order(w)
        for g in G:
            g.change_order(best_order)
        print("Chose order: " + str(w))

def spol(f, g):
    '''
    S-polynomial of f and g.
    '''
    R = f.parent()
    l = R.monomial_lcm(f.lm(), g.lm())
    return R(l / f.lt()) * f - R(l / g.lt()) * g

def buchberger(I, iterations=15, period=10):
    '''
    Very naive implementation of Buchberger's algorithm with support for a
    dynamic engine.
    '''
    G = [ BasisElement(g) for g in I.gens() ]
    P = [ (i, j) for i in range(len(G)) for j in range(i) ]
    dynamic = DynamicEngine(I.ring())
    while P:
        (i, j) = P[0]
        P = P[1:]
        s = spol(G[i].polynomial(), G[j].polynomial())
        f = s.reduce([ g.polynomial() for g in G])
        if f != 0:
            P = P + [ (i, len(G)) for i in range(len(G)) ]
            G.append(BasisElement(f))
            dynamic.next(G, iterations, period)
    J = ideal([ g.polynomial() for g in G ]).interreduced_basis()
    return len(J)

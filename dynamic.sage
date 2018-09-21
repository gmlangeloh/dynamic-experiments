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
            n = self.polynomial.parent().ngens()
            self._newton_polyhedron = self.polynomial().newton_polytope() + \
                negative_orthant(n)
        return self._newton_polyhedron

    def change_ring(self, R):
        self._polynomial = R(self._polynomial)

class DynamicEngine:

    def __init__(self, R):
        n = R.ngens()
        self._order = TermOrder("wdegrevlex", [1] * n)
        self._n = n
        self._ring = R

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
            vertex = [ lp.get_values(var)[i] for i in range(self._n) ]
            minkowski_vertex += vector(vertex)
            vertex_decomposition.append(vector(vertex))
        return w, minkowski_vertex, vertex_decomposition

    def _monomial_from_vertex(self, v):
        return prod([ R.gens()[i]^v[i] for i in self._n ])

    def next(self, G):
        w, v, v_decomposed = _random_minkowski_vertex(self, G)
        #Continue here
        #Unless I'm mistaken, we can also implement the restricted algorithm
        #easily working with Minkowski sums: it is enough to sum the previously
        #accepted vertex with the polytope of the new basis element!

def spol(f, g):
    '''
    S-polynomial of f and g.
    '''
    R = f.parent()
    l = R.monomial_lcm(f.lm(), g.lm())
    return R(l / f.lt()) * f - R(l / g.lt()) * g

def buchberger(I):
    '''
    Very naive implementation of Buchberger's algorithm with support for a
    dynamic engine.
    '''
    G = [ BasisElement(g) for g in I.gens() ]
    P = [ (i, j) for i in range(len(G)) for j in range(i) ]
    while P:
        (i, j) = P[0]
        P = P[1:]
        s = spol(G[i].polynomial(), G[j].polynomial())
        f = s.reduce([ g.polynomial() for g in G])
        if f != 0:
            P = P + [ (i, len(G)) for i in range(len(G)) ]
            G.append(BasisElement(f))
    return G

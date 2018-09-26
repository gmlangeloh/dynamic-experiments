'''
Simple implementation of Dynamic (unrestricted) Buchberger.
'''

from random import randint
from random import shuffle

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
        self._normal_fan = None
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

    def normal_fan(self):
        if self._normal_fan is None:
            fan = {}
            NP = self.newton_polyhedron()
            rays = [ -ieq.A() for ieq in NP.inequalities() ]
            for vertex in NP.vertices():
                indices = [ ieq.index() for ieq in vertex.incident() ]
                v_rays = [ rays[i] for i in indices ]
                fan[tuple(vertex)] = Cone(v_rays)
            self._normal_fan = fan
        return self._normal_fan

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
        self._call = -1

    def order(self):
        return _order

    def _random_vector(self):
        coords = [ randint(1, 1000) for i in xrange(self._n)]
        return vector(coords)

    def _find_vertices(self, G, w):
        #This is O(mk), where m = |G| and k = max(len(g)) for g in G
        #TODO definitely not efficient. Make this better
        vertices = []
        w_vec = vector(w)
        for g in G:
            NP = g.newton_polyhedron()
            min_val = float("inf")
            min_vertex = None
            for v in NP.vertices():
                val = vector(v) * w_vec
                if val < min_val:
                    min_val = val
                    min_vertex = v
            vertices.append(min_vertex)
        return vertices

    def _random_minkowski_vertex(self, G):
        '''
        Returns a random vertex of the Minkowski sum of the Newton polyhedra
        in G.
        '''
        w = self._random_vector()
        vertex_decomposition = self._find_vertices(G, w)
        v = sum([ vector(v) for v in vertex_decomposition])
        return list(w), tuple(v), vertex_decomposition

    def _monomial_from_vertex(self, v):
        return prod([ self._ring.gens()[i]^v[i] for i in xrange(self._n) ])

    def _ideal_from_decomposition(self, decomposition):
        G = [ self._monomial_from_vertex(v) for v in decomposition ]
        return ideal(G)

    def change_order(self, w, G):
        print("Chose order: " + str(w))
        self._order = TermOrder("wdegrevlex", w)
        for g in G:
            g.change_order(w)

    def next(self, G, iterations, period):
        '''
        Obtain a new monomial order for G.
        '''
        self._call += 1
        if self._call % period != 0:
            return
        for i in xrange(iterations):
            w, v, v_decomposed = self._random_minkowski_vertex(G)
            I = MonomialIdeal(self._ideal_from_decomposition(v_decomposed), w)
            if self._best_ideal is None or I < self._best_ideal:
                self._best_ideal = I
        best_order = self._best_ideal.weights()
        print("Chose order: " + str(best_order))
        self.change_order(best_order, G)

    def _neighborhood(self, G, v_decomposed):
        N = []
        for i in xrange(len(G)):
            NP = G[i].newton_polyhedron()
            NFan = G[i].normal_fan()
            v = v_decomposed[i]
            C1 = NFan[tuple(v)]
            for u in v.neighbors():
                if u.is_vertex():
                    C2 = NFan[tuple(u)]
                    C3 = C1.intersection(C2)
                    edge_normal = sum(C3)
                    if edge_normal == 0:
                        continue
                    new_decomposition = self._find_vertices(G, edge_normal)
                    new_decomposition[i] = u
                    new_vertex = sum([ vector(v) for v in new_decomposition ])
                    N.append((list(edge_normal), tuple(new_vertex), new_decomposition))
        return N

    def _local_search(self, G, iterations):
        visited = {}
        initial_sol = self._random_minkowski_vertex(G)
        i = 0
        to_visit = [ initial_sol ]
        current = None
        while i < iterations and to_visit:
            #Visit node
            w, v, v_decomposed = to_visit[0]
            to_visit = to_visit[1:]
            if v in visited:
                continue
            visited[v] = True
            I = MonomialIdeal(self._ideal_from_decomposition(v_decomposed), w)
            if current is None or I < current:
                current = I
            #Generate neighbors, put them in the queue
            to_visit += self._neighborhood(G, v_decomposed)
            i += 1
        if self._best_ideal is None or current < self._best_ideal:
            self._best_ideal = current
        best_order = self._best_ideal.weights()
        self.change_order(best_order, G)

    def _neighborhood2(self, G, v_decomposed):
        #TODO can do this MUCH better. Prioritize neighbors in new polys!
        for i in xrange(len(G) - 1, -1, -1):
            NP = G[i].newton_polyhedron()
            NFan = G[i].normal_fan()
            v = v_decomposed[i]
            C1 = NFan[tuple(v)]
            for u in v.neighbors():
                if u.is_vertex():
                    C2 = NFan[tuple(u)]
                    C3 = C1.intersection(C2)
                    edge_normal = sum(C3)
                    new_decomposition = self._find_vertices(G, edge_normal)
                    new_decomposition[i] = u
                    new_vertex = sum([ vector(v) for v in new_decomposition ])
                    yield (list(edge_normal), tuple(new_vertex), new_decomposition)

    def _local_search2(self, G, iterations):
        w, v, v_decomposed = self._random_minkowski_vertex(G)
        print(w)
        current = MonomialIdeal(self._ideal_from_decomposition(v_decomposed), w)
        to_visit = self._neighborhood2(G, v_decomposed)
        i = 0
        while i < iterations:
            try:
                #Visit node
                w, v, v_decomposed = next(to_visit)
                I = MonomialIdeal(self._ideal_from_decomposition(v_decomposed), w)
                if I < current:
                    print(w)
                    current = I
                    #Change current neighborhood to the next one (first improvement)
                    to_visit = self._neighborhood2(G, v_decomposed)
                i += 1
            except StopIteration:
                break
        #TODO make better checks if we should take current or not!!!
        self._best_ideal = current
        best_order = self._best_ideal.weights()
        self.change_order(best_order, G)

    def next2(self, G, iterations, period):
        self._call += 1
        if self._call % period != 0:
            return
        self._local_search(G, iterations)

    def next3(self, G, iterations, period):
        self._call += 1
        if self._call % period != 0:
            return
        self._local_search2(G, iterations)

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
    P = [ (i, j) for i in xrange(len(G)) for j in xrange(i) ]
    dynamic = DynamicEngine(I.ring())
    while P:
        (i, j) = P[0]
        P = P[1:]
        s = spol(G[i].polynomial(), G[j].polynomial())
        f = s.reduce([ g.polynomial() for g in G])
        if f != 0:
            P = P + [ (i, len(G)) for i in xrange(len(G)) ]
            G.append(BasisElement(f))
            dynamic.next3(G, iterations, period)
    J = ideal([ g.polynomial() for g in G ]).interreduced_basis()
    return len(J)

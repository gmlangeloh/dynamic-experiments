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
        self._graph = None
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
                fan[tuple(vertex)] = Cone(v_rays) #sum(v_rays)
            self._normal_fan = fan
        return self._normal_fan

    def graph(self):
        if self._graph is None:
            self._graph = self.newton_polyhedron().vertex_graph()
        return self._graph

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
    def __init__(self, I, w, vertices = None):
        self._ideal = I
        self._ring = I.ring()
        self._weights = w
        self._n = I.ring().ngens()
        self._vertices = vertices
        self._hilbert = None

    def hilbert(self):
        if self._hilbert is None:
            self._hilbert = self._ideal.hilbert_polynomial()
        return self._hilbert

    def vertices(self):
        return self._vertices

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

    def _lm_from_weights(self, g):
        max_val = None
        max_mon = None
        for v in g.graph().vertices():
            val = vector(v) * self.weights()
            if max_val is None or val > max_val:
                max_val = val
                max_mon = v
        return prod([ self._ring.gens()[i]^max_mon[i] \
                      for i in range(self._n) ]), max_mon

    def update(self, g):
        old_gens = self._ideal.gens()
        new_gen, new_v = self._lm_from_weights(g)
        self._ideal = ideal(old_gens + [ new_gen ])
        self._vertices += [ new_v ]
        self._hilbert = None

class DynamicEngine:

    def __init__(self, R):
        n = R.ngens()
        self._order = TermOrder("wdegrevlex", [1] * n)
        self._n = n
        self._ring = R
        self._best_ideal = None
        self._call = -1
        self._cone = None

    def order(self):
        return _order

    def _random_vector(self):
        coords = [ randint(1, 1000) for i in xrange(self._n)]
        return vector(coords)

    def _find_vertices(self, G, w):
        #This is O(mk), where m = |G| and k = max(len(g)) for g in G
        vertices = []
        w_vec = vector(w)
        for g in G:
            graph = g.graph()
            min_val = float("inf")
            min_vertex = None
            for v in graph.vertices():
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
        return list(w), vertex_decomposition

    def _monomial_from_vertex(self, v):
        return prod([ self._ring.gens()[i]^v[i] for i in xrange(self._n) ])

    def _ideal_from_decomposition(self, decomposition):
        G = [ self._monomial_from_vertex(v) for v in decomposition ]
        return ideal(G)

    def change_order(self, w, G):
        print(self._best_ideal.hilbert().degree(), self._best_ideal.hilbert().lc())
        self._order = TermOrder("wdegrevlex", w)
        for g in G:
            g.change_order(w)

    def _neighborhood(self, G, v_decomposed):
        for i in xrange(len(G) - 1, -1, -1):
            NFan = G[i].normal_fan()
            graph = G[i].graph()
            v = v_decomposed[i]
            for u in graph.neighbors(v):
                C = NFan[tuple(u)] #This could maybe be optimized...
                new_w = sum(C)
                new_decomposition = self._find_vertices(G, new_w)
                new_decomposition[i] = u
                yield (list(new_w), new_decomposition)

    def _local_search(self, G, iterations):
        w, v_decomposed = self._random_minkowski_vertex(G)
        current = MonomialIdeal(self._ideal_from_decomposition(v_decomposed),\
                                w, v_decomposed)
        to_visit = self._neighborhood(G, v_decomposed)
        i = 0
        while i < iterations:
            try:
                #Visit node
                w, v_decomposed = next(to_visit)
                I = MonomialIdeal(self._ideal_from_decomposition(v_decomposed),\
                                  w, v_decomposed)
                if I < current:
                    current = I
                    #Change current neighborhood to the next one (first improvement)
                    to_visit = self._neighborhood(G, v_decomposed)
                i += 1
            except StopIteration:
                break
        if self._best_ideal is None:
            self._best_ideal = current
            best_order = self._best_ideal.weights()
            self._cone = self._initial_cone(G, self._best_ideal.vertices())
            print("U: Chose order: " + str(w))
            self.change_order(best_order, G)
        else:
            #This can be optimized a little
            self._best_ideal.update(G[-1])
            if current < self._best_ideal:
                self._best_ideal = current
                self._cone = self._initial_cone(G, self._best_ideal.vertices())
                best_order = self._best_ideal.weights()
                print("U: Chose order: " + str(w))
                self.change_order(best_order, G)

    def _initial_cone(self, G, decomposition):
        C = G[0].normal_fan()[tuple(decomposition[0])]
        for i in xrange(1, len(G)):
            next_C = G[i].normal_fan()[tuple(decomposition[i])]
            C = C.intersection(next_C)
        return C

    def _recompute_basis_pairs(self, G):
        J = ideal([ g.polynomial() for g in G ]).interreduced_basis()
        G = [ BasisElement(g) for g in J ]
        P = [ (i, j) for i in xrange(len(G)) for j in xrange(i) ] #in unrestricted case, this is needed!
        return G, P

    def next(self, G, P, iterations, period):
        self._call += 1
        if self._call % period != 0:
            #Apply restricted algorithm most of the time!
            self._restricted_search(G)
            P = P + [ (i, len(G) - 1) for i in xrange(len(G)) ]
            return G, P
        else:
            #Sometimes, apply unrestricted search.
            self._local_search(G, iterations)
            return self._recompute_basis_pairs(G)

    def _restricted_search(self, G):
        '''
        A very naive implementation of Caboara and Perry's restricted
        dynamic step.
        '''
        NP = G[-1].newton_polyhedron()
        NFan = G[-1].normal_fan()
        min_cone = None
        current = None
        v_decomposition = self._best_ideal.vertices()
        for v in NP.vertices():
            C = NFan[tuple(v)]
            if self._cone is None:
                C_int = C
            else:
                C_int = C.intersection(self._cone)
            if not C_int.rays(): #If cone is empty
                continue
            w = sum(C_int)
            decomp = v_decomposition + [v]
            I = MonomialIdeal(self._ideal_from_decomposition\
                              (decomp), w, decomp)
            if current is None or I < current:
                current = I
                min_cone = C_int
        if current is None:
            print(G[-1].polynomial())
        self._cone = min_cone
        #print(self._best_ideal._ideal)
        #print(current._ideal)
        self._best_ideal = current
        best_order = self._best_ideal.weights()
        print("R: Chose order " + str(list(best_order)))
        self.change_order(list(best_order), G)

def spol(f, g):
    '''
    S-polynomial of f and g.
    '''
    R = f.parent()
    l = R.monomial_lcm(f.lm(), g.lm())
    return R(l / f.lt()) * f - R(l / g.lt()) * g

def buchberger(I, use_dynamic = True, iterations=15, period=10):
    '''
    Very naive implementation of Buchberger's algorithm with support for a
    dynamic engine.
    '''
    G = [ BasisElement(g) for g in I.gens() ]
    P = [ (i, j) for i in xrange(len(G)) for j in xrange(i) ]
    if use_dynamic:
        dynamic = DynamicEngine(I.ring())
    while P:
        (i, j) = P[0]
        P = P[1:]
        s = spol(G[i].polynomial(), G[j].polynomial())
        f = s.reduce([ g.polynomial() for g in G])
        if f != 0:
            #P = P + [ (i, len(G)) for i in xrange(len(G)) ]
            G.append(BasisElement(f))
            if use_dynamic:
                G, P = dynamic.next(G, P, iterations, period)
    J = ideal([ g.polynomial() for g in G ]).interreduced_basis()
    return len(J)

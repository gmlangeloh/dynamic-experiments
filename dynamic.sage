'''
Simple implementation of Dynamic (unrestricted) Buchberger.
'''

from random import randint
from random import shuffle

load("monomialideal.sage")

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
        coords = [ randint(1, 10) for i in xrange(self._n)]
        return vector(coords)

    def _not_so_random_vector(self):
        coords = [ randint(1, 10) for i in xrange(self._n)]
        perturb = randint(0, self._n - 1)
        coords[perturb] = randint(100, 10000)
        return vector(coords)

    def _find_vertices(self, G, w):
        #This is O(mk), where m = |G| and k = max(len(g)) for g in G
        vertices = []
        w_vec = vector(w)
        for g in G:
            graph = g.graph()
            max_val = float("-inf")
            max_vertex = None
            for v in graph.vertices():
                val = vector(v) * w_vec
                if val > max_val:
                    max_val = val
                    max_vertex = v
            vertices.append(max_vertex)
        return vertices

    def _random_minkowski_vertex(self, G):
        '''
        Returns a random vertex of the Minkowski sum of the Newton polyhedra
        in G.
        '''
        #w = self._random_vector()
        w = self._not_so_random_vector()
        print(w)
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
            g.change_order(list(w))

    def _neighborhood(self, G, v_decomposed):

        def diff(L1, L2):
            assert(len(L1) == len(L2))
            equal = 0
            for i in range(len(L1)):
                if L1[i] == L2[i]:
                    equal += 1
            return equal

        def mink(G):
            S = G[0].newton_polyhedron()
            for i in xrange(1, len(G)):
                S += G[i].newton_polyhedron()
            return S

        def vertex_in(v, P):
            m = len(v)
            for u in P.vertices():
                if all([v[i] == u[i] for i in range(m)]):
                    return u
            return None

        R = xrange(len(G) - 1, -1, -1)
        #M = mink(G)
        #R = xrange(len(G))
        for i in R:
            NFan = G[i].normal_fan()
            graph = G[i].graph()
            v = v_decomposed[i]
            #S = sum(vector(u) for u in v_decomposed)
            #S = sum([vector(u) for u in v_decomposed])
            #S = vertex_in(S, M)
            #s = len([u for u in S.neighbors() if u.is_vertex()])
            #print("N:" + str(s))
            Cv = NFan[tuple(v)]
            for u in graph.neighbors(v):
                Cu = NFan[tuple(u)] #This could maybe be optimized...
                C = Cv.intersection(Cu)
                new_w = 10000 * sum(C) + sum(Cu) #Integer version of taking an epsilon...
                new_decomposition = self._find_vertices(G, new_w)
                #S2 = sum(vector(u2) for u2 in new_decomposition)
                #S2 = sum([vector(u2) for u2 in new_decomposition])
                #S2 = vertex_in(S2, M)
                #if S2 in list(S.neighbors()):
                #    print("HEEY")
                assert(new_decomposition[i] == u)
                assert(Cu.relative_interior_contains(new_w))
                #e = diff(v_decomposed, new_decomposition)
                #print(e, len(v_decomposed) - e, sum([ vector(v) for v in new_decomposition]))
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
        #TODO for some reason, this always goes towards grevlex
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
        if any([ w == 0 for w in current.weights() ]):
            return
        self._cone = min_cone
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
    a = R.monomial_quotient(l, f.lt())
    b = R.monomial_quotient(l, g.lt())
    return a * f - b * g

reductions = 0
useful_reductions = 0

def clear_metadata():
    global reductions, useful_reductions
    reductions = 0
    useful_reductions = 0

def print_metadata(G):
    print("S-reductions: " + str(reductions))
    print("S-reductions (useful): " + str(useful_reductions))
    print("Size (before interreduction): " + str(len(G)))
    J = ideal([ g.polynomial() for g in G ]).interreduced_basis()
    print("Size (after interreduction): " + str(len(J)))

def buchberger(I, use_dynamic = True, iterations=15, period=10):
    '''
    Very naive implementation of Buchberger's algorithm with support for a
    dynamic engine.
    '''
    global reductions, useful_reductions
    clear_metadata()
    G = [ BasisElement(g) for g in I.gens() ]
    P = [ (i, j) for i in xrange(len(G)) for j in xrange(i) ]
    if use_dynamic:
        dynamic = DynamicEngine(I.ring())
    while P:
        (i, j) = P[0]
        P = P[1:]
        s = spol(G[i].polynomial(), G[j].polynomial())
        f = s.reduce([ g.polynomial() for g in G])
        reductions += 1
        if f != 0:
            useful_reductions += 1
            G.append(BasisElement(f))
            if use_dynamic:
                G, P = dynamic.next(G, P, iterations, period)
            else:
                P = P + [ (i, len(G) - 1) for i in xrange(len(G)) ]
    print_metadata(G)

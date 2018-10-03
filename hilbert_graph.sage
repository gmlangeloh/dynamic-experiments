r'''
Script to compute and experiment with the Hilbert graph of a polynomial
system.

The vertices of the Hilbert graph are those of the (affine) Newton polyhedron
of the system, defined as the Minkowski sum of the Newton polyhedra of its
polynomials, and the arcs are the neighborhood relations, oriented by the
decreasing direction of the Hilbert heuristic on the vertices.
'''

load("monomialideal.sage")

def negative_orthant(n):
    return Polyhedron(rays=(-identity_matrix(n)).rows())

def newton_polyhedron(p):
    n = p.parent().ngens()
    return p.newton_polytope() + negative_orthant(n)

def minkowski_sum(L):
    S = newton_polyhedron(L[0])
    for i in xrange(1, len(L)):
        S += newton_polyhedron(L[i])
    return S

def normal_cone(V, P):
    rays = [ -ieq.A() for ieq in P.inequalities() ]
    indices = [ ieq.index() for ieq in V.incident() ]
    v_rays = [ rays[i] for i in indices ]
    return Cone(v_rays)

def ideal_lt_from_vertex(V, P, G):
    C = normal_cone(V, P)
    normal = vector(sum(C))
    LMs = []
    for g in G:
        max_val = -float("inf")
        lm = None
        for m in g.monomials():
            val = normal * vector(m.exponents()[0])
            if val > max_val:
                max_val = val
                lm = m
        LMs.append(lm)
    I = ideal(LMs)
    return MonomialIdeal(I)

#TODO be careful when using sinks... there may exist sink-multivertices
class HilbertGraph:

    def __init__(self, polyhedron, G):
        gr = polyhedron.vertex_graph()
        self._polyhedron = polyhedron
        self._polynomial_system = G
        self._graph = DiGraph()
        self._graph.add_vertices(gr.vertices())
        edges = []
        for v in gr.vertices():
            Iv = ideal_lt_from_vertex(v, polyhedron, G)
            for u in v.neighbors():
                if u.is_vertex():
                    Iu = ideal_lt_from_vertex(u, polyhedron, G)
                    if Iv < Iu:
                        e = [u, v]
                        edges.append(e)
                    elif Iv > Iu:
                        e = [v, u]
                        edges.append(e)
                    else:
                        e1 = [u, v]
                        e2 = [v, u]
                        edges.append(e1)
                        edges.append(e2)
        self._graph.add_edges(edges)

    def graph(self):
        return self._graph

    def heights(self):

        def height(v):
            Iv = ideal_lt_from_vertex(v, self._polyhedron, self._polynomial_system)
            return Iv.hilbert().degree() * 10 + float(Iv.hilbert().lc())

        h_dict = {}
        for v in self.graph().vertices():
            h = height(v)
            if h not in h_dict:
                h_dict[h] = [v]
            else:
                h_dict[h].append(v)
        return h_dict

def graph(G):
    r'''
    Computes the Hilbert graph of the polynomial system G.
    '''
    P = minkowski_sum(G)
    return HilbertGraph(P, G)

def graph_partial(G, i):
    P = minkowski_sum([G[i]])
    return HilbertGraph(P, [G[i]])

def test():
    R = PolynomialRing(GF(32003), "x", 3)
    G = [ R.random_element(degree=5, terms=10) for j in range(3) ]
    gr = graph(G)
    gr._graph.show(edge_style="dashed")

'''
Useful module to test the effects of making an initial analysis

We can print:
- best/worst Hilbert Functions, average
- compare with the ordering Perry's algorithm would choose in the same number of steps
- compare with grevlex as well
'''

def minkowski(I):
    n = I.ring().ngens()
    result = Polyhedron(rays=(-identity_matrix(n)).rows())
    for g in I.gens():
        result += g.newton_polytope()
    return result

def normal_cone(v, P):
    inequalities = P.inequalities()
    indices = [ ieq.index() for ieq in v.incident() ]
    rays = [ -inequalities[i].A() for i in indices ]
    return Cone(rays)

def normal(v, P):
    return sum(normal_cone(v, P))

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
    HP = ideal(LMs).hilbert_polynomial()
    return G, w, HP.degree(), HP.lc()

def evaluate_all(I):
    Mink = minkowski(I)
    mindeg = float("inf")
    maxdeg = float("-inf")
    for v in Mink.vertex_generator():
        G, w, deg, coef = evaluate(v, Mink, I)
        LMs = [ g.lm() for g in G]
        #print LMs, w, deg, coef
        if deg < mindeg:
            mindeg = deg
        if deg > maxdeg:
            maxdeg = deg
    print mindeg, maxdeg

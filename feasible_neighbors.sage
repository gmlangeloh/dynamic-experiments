'''
Check if neighbors in Minkowski sum are feasible with decent probability or not.
'''

from random import randint, seed
import fcntl
import sys

load("benchmarks.sage")

def feasible(summands, all_vertices):
    '''
    Returns true iff sum(summands) is a vertex of the Minkowski sum of
    the polyhedra with vertices given by `all_vertices`

    - summands -- a list of vertices, one of each summand
    - all_vertices -- a list of lists of vertices, each list of vertices
    corresponding to a distinct polyhedron
    '''

    #Build linear program from summands
    n = len(summands[0])
    k = len(all_vertices)
    from sage.numerical.backends.generic_backend import get_solver
    import sage.numerical.backends.glpk_backend as backend
    lp = get_solver(solver="GLPK")
    lp.solver_parameter("simplex_or_intopt", "simplex_only")
    E = 0.001
    lp.add_variables(n, lower_bound=E)

    for i in xrange(k):
        V = summands[i]
        vertices = all_vertices[i]
        for W in vertices:
            if V == W:
                continue
            coefs = list(vector(V) - vector(W))
            lp.add_linear_constraint(list(zip(xrange(n), coefs)), E, None)

    lp.set_objective([1.0] * n)
    lp.set_sense(-1)

    try:
        lp.solve()
    except Exception as e:
        s = str(e)
        if "no feasible" in s:
            return False
        else:
            raise e

    return True

def has_normal(vertices, w):
    '''
    Returns the vertex from `vertices` that has `w` as normal vector.
    Assumes one of them does (this work in the case of Newton polyhedra
    as long as w > 0, which is what I want to test anyway)
    '''

    w = vector(w)
    max_val = 0.0
    max_vertex = None
    for v in vertices:
        val = w * vector(v)
        if val > max_val:
            max_val = val
            max_vertex = v
    return max_vertex

def initial_normal(n, random = True):
    '''
    Returns a initial normal vector that can be random or just the
    [1, 1, ..., 1] vector.
    '''

    w = [1] * n
    if random:
        w = [ randint(1, 60000) for i in xrange(n) ]
    return w

def negative_orthant(n):
    return Polyhedron(rays=(-identity_matrix(n)).rows())

def count_neighbors(G):

    #Build Newton polyhedra
    n = G[0].parent().ngens()
    k = len(G)
    NPs = [ g.newton_polytope() + negative_orthant(n) for g in G ]
    all_vertices = [ NP.vertices() for NP in NPs ]

    #Find initial vertex in Minkowski sum
    w = initial_normal(n)
    init_summands = [ has_normal(all_vertices[i], w) for i in xrange(k) ]
    summands = list(init_summands)
    assert feasible(summands, all_vertices), "Infeasible initial solution... bug?"

    #Count neighbors of the initial vertex in Minkowski sum, compare to amount
    #of potential neighbors
    num_neighbors = 0
    num_total = 0
    for idx, polyhedron in enumerate(all_vertices):
        for vertex in polyhedron:
            summands[idx] = vertex
            if feasible(summands, all_vertices):
                num_neighbors += 1
            num_total += 1
        summands = list(init_summands)

    return num_neighbors, num_total

def neighbor_data(instance):

    G = Benchmark(instance).ideal.gens()
    name = instance.split("/")[-1].split(".")[0]
    if G[0].parent().ngens() > 8:
        return

    sys.stderr.write("Starting: " + name + "\n")

    seed(0)
    for i in xrange(1, 31): #Number of repetitions
        neighbors, total = count_neighbors(G)
        fcntl.lockf(sys.stdout, fcntl.LOCK_EX) #only one process writes at a time
        print name, G[0].parent().ngens(), neighbors, total
        fcntl.lockf(sys.stdout, fcntl.LOCK_UN)

    sys.stderr.write("Finished: " + name + "\n")

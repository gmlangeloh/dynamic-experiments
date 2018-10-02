r'''
Experiments involving the complexity of common refinements compared to
Minkowski sums.

Result: common refinement, in Sage, is MORE expensive than Minkowski
sums. This may be a consequence of the lattice structure that Sage
automatically computes for the former that it only computes when needed
for the latter.
'''

from time import time

R = PolynomialRing(GF(32003), "x", 5)

def common_refinement(L):
    r'''
    The common refinement of a list of polynomials
    '''
    C = L[0].newton_polytope().normal_fan()
    for i in xrange(1, len(L)):
        C = C.common_refinement(L[i].newton_polytope().normal_fan())
    return C

def minkowski_sum(L):
    r'''
    The Minkowski sum of a list of polynomials
    '''
    S = L[0].newton_polytope()
    for i in xrange(1, len(L)):
        S += L[i].newton_polytope()
    return S

for i in range(1, 20):
    L = [ R.random_element(terms=50) for j in range(i) ]
    t_ref = time()
    C = common_refinement(L)
    t_ref = time() - t_ref
    t_sum = time()
    S = minkowski_sum(L)
    t_sum = time() - t_sum
    print("Testing " + str(i) + " polytopes...")
    print("Common refinement: " + str(t_ref))
    print("Minkowski sum: " + str(t_sum))

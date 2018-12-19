# cython: profile = False
# cython: boundscheck = False
# cython: wraparound = False
# clang c++
# cinclude $SAGE_ROOT/local/include/singular
# clib m readline Singular givaro gmpxx gmp

r"""
    This Sage/Cython code implements a dynamic algorithm to compute a Groebner basis,
    based on the Dynamic Buchberger Algorithm of Caboara.
    It uses a comparison with boundary vectors to minimize the number of polynomials
    tested.
    This is merely a study implementation.

    USAGE:

        see the function dynamic_gb
"""

from types cimport *

from copy import copy
from random import randint, choice

import cython
from cython.parallel cimport prange

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix

from sage.misc.misc_c import prod

from sage.rings.infinity import Infinity

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder

from sage.rings.integer_ring import IntegerRing

from sage.modules.free_module_element import vector
from sage.rings.real_double import RDF
from sage.geometry.polyhedron.constructor import Polyhedron

from sage.functions.other import floor, ceil

# record keeping
cdef int rejections, monomials_eliminated, number_of_programs_created

cpdef GLPKBackend make_solver(int n):
  r"""
  Creates an empty model in a GLPK backend solver.

  OUTPUTS:
  - a GLPK backend
  """

  from sage.numerical.backends.generic_backend import get_solver
  import sage.numerical.backends.glpk_backend as backend
  lp = get_solver(solver="GLPK")
  lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
  if n > 0:
    lp.add_variables(n)

  return lp

cpdef void append_linear_program(GLPKBackend glpk, MPolynomial_libsingular p, int k):
  r"""
  Appends constraints and variables of the Newton polyhedron of `p` to the current
  linear program in `glpk`.

  INPUTS:

  - `glpk` -- current representation of the linear programming model in GLPK
  - `p` -- polynomial whose affine Newton Polyhedron will be added to the linear programming model
  - `k` -- number of polynomials in current basis (including p)
  """

  cdef int n = p.parent().ngens()
  cdef tuple a
  cdef list variables, coefs

  #Remove previous auxiliary constraints
  if glpk.nrows() > 0:
    for i in xrange(glpk.nrows() - 1, glpk.nrows() - n - 1, -1):
      glpk.remove_constraint(glpk.nrows() - 1)

  NP = p.newton_polytope() + Polyhedron(rays=(-identity_matrix(n)).rows())
  lp = NP.to_linear_program(solver="GLPK")
  cdef int cols = glpk.ncols()

  glpk.add_variables(n)

  for lb, a, ub in lp.constraints():
    variables, coefs = a
    variables = [ i + cols - n for i in variables ]
    glpk.add_linear_constraint(list(zip(variables, coefs)), lb, ub)

  #Add auxiliary constraints
  for i in xrange(n):
    L = [ (i + j * n, -1) for j in xrange(k) ]
    L.append((i + cols, 1))
    glpk.add_linear_constraint(L, 0.0, 0.0)

  return

@cython.profile(True)
cpdef void init_linear_program(GLPKBackend lp, int n):
  lp.add_variables(n) #Add variables representing the negative orthant summand
  return

@cython.profile(True)
cpdef void update_linear_program(GLPKBackend lp, MPolynomial_libsingular p, list vertices):

  cdef int n = p.parent().ngens()
  NP = p.newton_polytope() + Polyhedron(rays=(-identity_matrix(n)).rows())
  Vs = NP.vertices()
  vertices.append((lp.ncols(), Vs))
  cdef int l = len(Vs)
  lp.add_variables(l)

  #Remove previous constraints determining the point of the Minkowski sum
  cdef list rows = []
  first = False
  if lp.nrows() >= n:
    for i in xrange(lp.nrows() - n, lp.nrows()):
      rows.append(lp.row(i))
  if lp.nrows() > n:
    lp.remove_constraints(xrange(lp.nrows() - n, lp.nrows()))
  else:
    first = True

  #Insert auxiliary constraint
  lp.add_linear_constraint(list(zip(xrange(lp.ncols() - l, lp.ncols()), [1]*l)), 1.0, 1.0)

  #Insert constraint determining the point of the Minkowski sum
  cdef list var_indices, var_coefs
  for i in xrange(n):
    if first:
      var_indices = [i, i + n] + range(lp.ncols() - l, lp.ncols())
      var_coefs = [-1.0, -1.0] + [ Vs[j][i] for j in xrange(l) ]
    else:
      var_indices = rows[i][0] + range(lp.ncols() - l, lp.ncols())
      var_coefs = rows[i][1] + [ Vs[j][i] for j in xrange(l) ]
    lp.add_linear_constraint(list(zip(var_indices, var_coefs)), 0.0, 0.0)

  return

@cython.profile(True)
cpdef list weight_vector(GLPKBackend lp, int n):
  r"""
  Returns the weight vector currently used as objective function in the linear
  programming model lp.

  INPUTS:

  - `lp` -- a GLPK linear programming model
  - `n` -- number of variables in the polynomial system

  OUTPUTS:

  - a list representing a weight vector
  """
  cdef list w = []
  cdef int i
  #for i in xrange(lp.ncols() - n, lp.ncols()):
  #  w.append(lp.objective_coefficient(i))
  for i in xrange(0, n):
    w.append(lp.objective_coefficient(i))
  return w

cpdef Vector_real_double_dense tableau_row(GLPKBackend lp, int i, list nonbasic):
  r"""
  Returns the i-th row of the simplex tableau of `lp` as a dense vector.

  INPUTS:

  - `lp` -- a previously solved GLPK linear programming model
  - `i` -- index of the desired row
  - `nonbasic` -- list of nonbasic variable indices of `lp`

  OUTPUT:

  - a vector representation of the i-th row of `lp`
  """

  cdef list row = [], row_indices, row_coeffs
  row_indices, row_coeffs = lp.eval_tab_row(i)
  cdef int idx = 0, m = lp.nrows(), j
  for j in nonbasic:
    if idx < len(row_indices) and j == row_indices[idx]:
      if j >= m:
        row.append(-row_coeffs[idx])
      else:
        row.append(row_coeffs[idx])
      idx += 1
    else:
      row.append(0.0)
  return vector(RDF, row)

@cython.profile(True)
cpdef list apply_sensitivity_range(float lower, float upper, GLPKBackend lp, int change_idx, int n):

  cdef float epsilon = 0.001
  if abs(lower - round(lower)) < epsilon:
    lower = float(round(lower))
  if abs(upper - round(upper)) < epsilon:
    upper = float(round(upper))

  cdef float increment
  if lower == float("-inf") and upper == float("inf"):
    return weight_vector(lp, n)
  if upper == float("inf"):
    increment = ceil(lower) - 1
    if lp.objective_coefficient(change_idx) + increment <= 0:
      return weight_vector(lp, n)
  elif lower == float("-inf"):
    increment = floor(upper) + 1
  else:
    #Lower coefficient, when possible, 50% of the time
    increment = ceil(lower) - 1
    if lp.objective_coefficient(change_idx) + increment <= 0:
      increment = floor(upper) + 1
    elif randint(0, 1):
      increment = floor(upper) + 1
      #pass

  assert abs(increment) > 0
  cdef float old_value = lp.objective_coefficient(change_idx)
  lp.objective_coefficient(change_idx, old_value + increment)

  return weight_vector(lp, n)

@cython.profile(True)
cpdef tuple sensitivity(GLPKBackend lp, int n, int k):
  r"""
  Changes the current lp objective function to point to a neighbor.

  INPUTS:

  - `lp` - a GLPK linear programming model pointing to the current chosen order
  - `n` - the number of variables in the input polynomial system
  - `k` - the number of polynomials currently in the basis
  """

  cdef list basic = [], nonbasic = []
  cdef int m = lp.nrows(), idx, i, j

  #Classify variables in basic/nonbasic
  for i in xrange(m):
    if lp.is_slack_variable_basic(i):
      basic.append(i)
    else:
      nonbasic.append(i)
  for i in xrange(lp.ncols()):
    if lp.is_variable_basic(i):
      basic.append(m + i)
    else:
      nonbasic.append(m + i)

  #Compute zN and DzN
  cdef int change = randint(0, n-1) #Choose which coefficient should be changed
  cdef list zN = []
  cdef Vector_real_double_dense DcN = vector(RDF, [ 0.0 ] * len(nonbasic))
  for idx, i in enumerate(nonbasic):
    if i >= m:
      zN.append(-lp.get_col_dual(i - m))
      if (i - m) % n == change:
        DcN[idx] = 1.0
    else:
      zN.append(lp.get_row_dual(i))

  cdef Vector_real_double_dense DzN
  if lp.is_variable_basic(change + k * n):
    DzN = tableau_row(lp, m + change + k * n, nonbasic) - DcN
  else:
    DzN = -DcN

  #Compute min/max interval
  cdef float t
  cdef float upper = float("inf")
  cdef float lower = float("-inf")

  cdef float epsilon = 0.00001
  for i in xrange(len(nonbasic)):
    if abs(DzN[i]) < epsilon and abs(zN[i]) < epsilon:
      t = 0.0
    elif abs(zN[i]) < epsilon:
      t = float("-inf") if DzN[i] > 0 else float("inf")
    else:
      t = -zN[i] / DzN[i]
    if t >= 0 and t < upper:
      upper = t
    if t <= 0 and t > lower:
      lower = t

  print "sensitivity:", lower, upper
  assert upper >= 0 and lower <= 0, "Inconsistent sensitivity range"

  #return apply_sensitivity_range(lower, upper, lp, n * k + change, n), n * k + change
  return [], n * k + change

cpdef void sensitivity_warm_start(GLPKBackend lp, GLPKBackend aux):

  cdef int i
  for i in xrange(aux.ncols() - 1):
    if lp.get_row_stat(i) != 1:
      aux.set_col_stat(i, 1)
    else:
      aux.set_col_stat(i, 2)
  for i in xrange(aux.nrows() - 1):
    if lp.get_col_stat(i) != 1:
      aux.set_row_stat(i, 1)
    else:
      aux.set_row_stat(i, 2)
  aux.set_row_stat(aux.nrows() - 1, 1)
  aux.set_col_stat(aux.ncols() - 1, 2)
  aux.warm_up()

cpdef list lp_column(GLPKBackend lp, int col_idx):

  cdef int i, j, idx
  cdef list column = [], row_indices, row_coeffs
  for i in xrange(lp.nrows()):
    row_indices, row_coeffs = lp.row(i)
    for idx, j in enumerate(row_indices):
      if j == col_idx:
        column.append((i, row_coeffs[idx]))
        break
  return column

cpdef list lp_bounds(GLPKBackend lp):

  cdef int i
  cdef list b = []
  for i in xrange(lp.nrows()):
    lb, ub = lp.row_bounds(i)
    if ub is not None:
      b.append(ub)
    elif lb is not None:
      b.append(lb)

  assert len(b) == lp.nrows()

  return b

@cython.profile(True)
cpdef list wide_sensitivity(GLPKBackend lp, int n, int coef = 0):
  r"""
  Implements the sensitivity analysis idea from Jensen et al, 1997.
  """

  #Make model here
  cdef int coef_change_idx, idx
  cdef int i
  if coef == 0:
    #coef_change_idx = randint(lp.ncols() - n, lp.ncols() - 1)
    #coef_change_idx = randint(0, n-1)
    m = 100000.0
    for i in xrange(n):
      if lp.objective_coefficient(i) < m:
        m = lp.objective_coefficient(i)
        idx = i
    coef_change_idx = idx
  else:
    coef_change_idx = coef
  cdef GLPKBackend glpk = make_solver(0)
  cdef int num_vars = lp.nrows()
  glpk.add_variables(num_vars, lower_bound=None) #Variables of the dual problem
  glpk.add_variable(lower_bound=None, upper_bound=None, binary=False, continuous=True,integer=False, obj=1.0) #The gamma variable
  cdef int gamma = glpk.ncols() - 1

  for i in xrange(lp.ncols()):
    if i != coef_change_idx:
      glpk.add_linear_constraint(lp_column(lp, i), lp.objective_coefficient(i), None)
    else:
      glpk.add_linear_constraint(lp_column(lp, i) + [(gamma, -1.0)], lp.objective_coefficient(i), None)

  cdef float epsilon = 0.0001
  cdef float zeta = lp.get_objective_value()
  #assert abs(zeta - c*x) < epsilon, "zeta is weird " + str(zeta) + " " + str(c*x)
  glpk.add_linear_constraint(list(zip(xrange(num_vars), lp_bounds(lp))) + [(gamma, -lp.get_variable_value(coef_change_idx))], zeta, zeta)

  #Solve for maximization
  cdef float upper
  glpk.set_sense(+1)
  #lp.write_lp("wtf.lp")
  #glpk.write_lp("sensitivity.lp")
  sensitivity_warm_start(lp, glpk)
  try:
    glpk.solve()
    upper = glpk.get_variable_value(gamma)
  except Exception as e:
    s = str(e)
    if "no feasible" in s:
      upper = 0.0
    elif "unbounded" in s:
      upper = float("inf")
    else:
      raise e

  #Solve for minimization
  cdef float lower
  glpk.set_sense(-1)
  sensitivity_warm_start(lp, glpk)
  try:
    glpk.solve()
    lower = glpk.get_variable_value(gamma)
  except Exception as e:
    s = str(e)
    if "no feasible" in s:
      lower = 0.0
    elif "unbounded" in s:
      lower = float("-inf")
    else:
      raise e

  print "wide sensitivity:", lower, upper
  assert lower < epsilon and upper > -epsilon, "Inconsistent sensitivity range"

  return apply_sensitivity_range(lower, upper, lp, coef_change_idx, n)

@cython.profile(True)
cpdef list find_monomials2(GLPKBackend lp, MPolynomialRing_libsingular R, list vertices, int k):
  #vertices is a list with tuples (idx, tuple) where tuple is a tuple with vertices, and idx is the index of the
  #first variable in lp referring to these vertices
  cdef int n = R.ngens()
  cdef int l, i, j
  cdef float val, epsilon = 0.0001
  cdef list LTs = []
  cdef MPolynomial_libsingular LM
  for i in xrange(k):
    idx = vertices[i][0]
    l = len(vertices[i][1])
    for j in xrange(l):
      val = lp.get_variable_value(idx + j)
      if abs(val - 1.0) < epsilon:
        break
    vertex = vertices[i][1][j]
    LM = prod([ R.gens()[j]**vertex[j] for j in xrange(n)])
    LTs.append(LM)

  return LTs

@cython.profile(True)
cpdef list find_monomials(GLPKBackend lp, MPolynomialRing_libsingular R, int k):
  r"""
  Obtains the leading monomials chosen by the order w in the linear programming model.

  INPUTS:

  - `lp` -- linear programming model (must be already solved)
  - `R` -- a polynomial ring
  - `k` -- number of polynomials in current basis

  OUTPUTS:

  - a list of leading monomials
  """
  cdef int n = R.ngens()
  cdef list LTs = []
  cdef int i, j, e
  cdef MPolynomial_libsingular monomial
  for i in xrange(k):
    #build i-th leading monomial
    monomial = R(1)
    for j in xrange(n):
      e = round(lp.get_variable_value(j + i * n))
      monomial *= R.gens()[j]**e
    LTs.append(monomial)
  return LTs

first = True
#TODO why does the number of iterations affect performance so much?
@cython.profile(True)
cpdef tuple choose_simplex_ordering(list G, list current_ordering, GLPKBackend lp, list vertices, int iterations = 5):
  r"""

  INPUTS:

  - `G` -- the current system of generators
  - `current_ordering` -- the current weight ordering
  - `lp` -- previous GLPK linear programming model
  - `vertices` -- TODO
  - `iterations` -- number of neighbors to visit

  OUTPUTS:

  - a list of weights representing a monomial order
  """
  global first
  cdef MPolynomialRing_libsingular R = G[0].value().parent()
  cdef MPolynomialRing_libsingular newR
  cdef int k = len(G)
  cdef int n = R.ngens()
  cdef int i, j, it = 0
  cdef list CLTs, LTs, oldLTs, w, best_w

  #TODO use iteration_limit, if we can keep the right value of w
  #lp.solver_parameter("iteration_limit", 2**31 - 1)

  #Initial random ordering
  if first:
    #w = [ randint(1, 10000) for i in xrange(n) ]
    init_linear_program(lp, n)
    #w = [10000.0] * n
    w = [4627, 8716, 1234, 876, 1038]
    first = False
  else:
    w = current_ordering
  best_w = w
  print w

  #Transform last element of G to linear program, set objective function given by w and solve
  #append_linear_program(lp, G[len(G)-1].value(), k)
  update_linear_program(lp, G[k-1].value(), vertices)
  #lp.set_objective([0] * n * k + w)
  lp.set_objective(w)
  lp.solve()

  #Get current LTs to compare with Hilbert heuristic
  newR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(w))
  LTs = find_monomials2(lp, newR, vertices, k)
  oldLTs = LTs
  CLTs = [ (newR.ideal(LTs).hilbert_polynomial(), newR.ideal(LTs).hilbert_series(), w ) ]

  #lp.solver_parameter("iteration_limit", 1)
  #Do sensitivity analysis to get neighbor, compare
  while it < iterations:
    #w, c = sensitivity(lp, n, k)
    w = wide_sensitivity(lp, n)
    lp.solve()
    newR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(w))
    LTs = find_monomials2(lp, newR, vertices, k)
    print [LTs[i] == oldLTs[i] for i in xrange(len(LTs))].count(False), len(LTs)
    if [LTs[i] == oldLTs[i] for i in xrange(len(LTs))].count(False) == 0:
      continue
    #elif lp.get_objective_value() == lp.best_known_objective_bound():
    #  print "before", w
    #  w = find_objective_dual(lp, n)
    #  print "after", w
    CLTs.append((newR.ideal(LTs).hilbert_polynomial(), newR.ideal(LTs).hilbert_series(), w))
    CLTs.sort(cmp=hs_heuristic)
    print "candidate 1:", CLTs[0][2], CLTs[0][0].degree(), CLTs[0][0].lc()
    print "candidate 2:", CLTs[1][2], CLTs[1][0].degree(), CLTs[1][0].lc()
    best_w = CLTs[0][2] #Take first improvement
    if best_w == w:
      oldLTs = LTs
    #lp.set_objective([0] * n * k + best_w)
    lp.set_objective(best_w)
    lp.solve()
    CLTs = CLTs[:1]
    it += 1

  ##Compare with current_ordering - keep the current one if they tie!
  #newR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(current_ordering))
  #LTs = [ newR(G[k].value()).lm() for k in xrange(len(G)) ]
  #CLTs.insert(0, (newR.ideal(LTs).hilbert_polynomial(), newR.ideal(LTs).hilbert_series(), current_ordering))
  #CLTs.sort(cmp=hs_heuristic)
  #best_w = CLTs[0][2]
  #lp.set_objective(best_w * k)
  #lp.solve()

  return best_w, vertices

@cython.profile(True)
cpdef list choose_random_ordering(list G, list current_ordering, int iterations = 10):
  r"""
  Chooses a weight vector for a term ordering for the basis ``G`` that is optimal
  with respect to the Hilbert tentative function on G among randomly generated orders.

  INPUTS:

  - `G` -- a basis of a polynomial ideal
  - `current_ordering` -- the current ordering of G, as a list of weights

  OUTPUTS:

  - a weighted ordering, as a list of weights
  """

  cdef int n = G[0].value().parent().ngens()
  cdef list rand_weights = [ current_ordering ]
  cdef list w, CLTs, LTs, best_order
  cdef MPolynomialRing_libsingular R = G[0].value().parent()
  cdef MPolynomialRing_libsingular newR

  cdef int i
  for i in xrange(iterations):
    #Choose random vector
    w = [ randint(1, 10) for i in xrange(n) ]
    rand_weights.append(w)

  #Compute CLTs
  CLTs = []
  for w in rand_weights:
    newR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(w))
    LTs = [ newR(G[i].value()).lm() for i in xrange(len(G)) ]
    CLTs.append((w, LTs))

  #Evaluate CLTs with Hilbert function
  best_order = min_weights_by_Hilbert_heuristic(R, CLTs)

  return best_order

@cython.profile(True)
cpdef list choose_local_ordering(list G, list current_ordering, int iterations = 50):
  r"""
  Chooses a weight vector for polynomial system `G` randomly and then optimizes it
  locally for a few iterations using small perturbations.

  INPUTS:

  - `G` -- a basis of a polynomial ideal
  - `current_ordering` -- the current ordering of G, as a list of weights
  - `iterations` -- number of neighbors to evaluate

  OUTPUTS:

  - a weighted ordering, as a list of weights
  """

  cdef int n = G[0].value().parent().ngens()
  cdef list curr_w, w, LTs, CLTs
  cdef int i, j, k, incr
  cdef MPolynomialRing_libsingular R = G[0].value().parent()
  cdef MPolynomialRing_libsingular newR

  #Choose random initial vector
  curr_w = [ randint(1, 1000) for i in xrange(n) ]
  print "initial order: ", curr_w
  newR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(curr_w))
  LTs = [ newR(G[k].value()).lm() for k in xrange(len(G)) ]
  CLTs = [ (newR.ideal(LTs).hilbert_polynomial(), newR.ideal(LTs).hilbert_series(), curr_w) ]

  #Compute perturbations
  for i in xrange(iterations):
    w = curr_w[:]
    j = randint(0, n-1)
    if w[j] < 2:
      incr = 1
    else:
      incr = choice([1, -1])
    w[j] += incr

    #Find LTs w.r.t current order w
    newR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(w))
    LTs = [ newR(G[k].value()).lm() for k in xrange(len(G)) ]
    CLTs.append((newR.ideal(LTs).hilbert_polynomial(), newR.ideal(LTs).hilbert_series(), w))

    #Choose best one among current and perturbed orders
    CLTs.sort(cmp=hs_heuristic)
    curr_w = CLTs[0][2] #Work in a first improvement basis
    CLTs = CLTs[:1]

  #Choose best order between `current_ordering` and `curr_w`
  newR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(current_ordering))
  LTs = [ newR(G[k].value()).lm() for k in xrange(len(G)) ]
  CLTs.append((newR.ideal(LTs).hilbert_polynomial(), newR.ideal(LTs).hilbert_series(), current_ordering))
  CLTs.sort(cmp=hs_heuristic)
  curr_w = CLTs[0][2]
  print "finally chose order: ", curr_w
  return curr_w

@cython.profile(True)
cpdef tuple choose_ordering_unrestricted(list G, list old_vertices):
  r"""
  Chooses a weight vector for a term ordering for the basis ``G`` that minimizes the Hilbert
  tentative function among the viable orderings from G. See [Gritzmann & Sturmfels 1993] for
  a description of the algorithm.

  INPUTS:

  - ``G`` -- a basis of a polynomial ideal
  - ``old_vertices`` -- list of tuples (``v``, summands) where ``v`` is a
  vertex of the previous Minkowski sum

  OUTPUTS:

  - a weighted ordering optimizing the Hilbert heuristic
  - a list of tuples (``v``, summands) where ``v`` is a vertex of the
  Minkowski sum
  """

  cdef list CLTs, new_vertices, lold, summands, gens, w, LTs
  cdef Vector_integer_dense vec_v1, vec_v2
  cdef tuple tup, constraint
  cdef MPolynomialRing_libsingular R = G[0].value().parent() # current ring
  cdef MPolynomial_libsingular p = G[len(G)-1].value()
  cdef int n = R.ngens()

  #STEP 1: find candidate LTs computing a Minkowski sum

  #Affine Newton polyhedron of the newest basis element
  new_polyhedron = p.newton_polytope() + Polyhedron(rays=(-identity_matrix(n)).rows())

  #List of tuples potentially in the new Minkowski sum
  #(vertex, summands)
  new_vertices = []

  #Compute the Minkowski sum
  for v1, lold in old_vertices:
    for v2 in new_polyhedron.vertex_generator():
      vec_v1 = vector(v1)
      vec_v2 = v2()
      new_vertices.append((list(vec_v1 + vec_v2), lold + [list(v2)]))

  #If this is the first time we are calling this, old_vertices is empty.
  if not old_vertices:
    new_vertices = [ (list(v),[list(v)]) for v in new_polyhedron.vertices() ]

  M = new_polyhedron.parent().element_class(new_polyhedron.parent(), \
                                            [[tup[0] for tup in new_vertices], new_polyhedron.rays(), []], \
                                            None)

  #Keep only vertices that are extreme points in the Minkowski sum ``M``
  new_vertices = [ tup for tup in new_vertices \
                   if tup[0] in map(list, M.vertices()) ]

  #STEP 2: evaluate candidates using the Hilbert heuristic

  #Compute the list of candidate LTs from the list of vertices
  CLTs = []
  for tup in new_vertices:
    summands = tup[1]
    gens = []
    for s in summands:
      gens.append(prod([ R.gens()[i]**s[i] for i in xrange(n) ]))
    CLTs.append(gens)
  LTs = min_CLT_by_Hilbert_heuristic(R, CLTs)

  #TODO can probably do this more efficiently with the normal cone model?
  #STEP 3: obtain a weight vector for the chosen order using linear programming

  import sage.numerical.backends.glpk_backend as glpk_backend
  lp = new_linear_program()
  lp.solver_parameter(glpk_backend.glp_simplex_or_intopt, glpk_backend.glp_simplex_then_intopt)

  # need positive weights
  for k in xrange(n):
    lp.add_constraint(lp[k],min=tolerance_cone)
    #lp.set_min(lp[k],tolerance_cone)
    lp.set_integer(lp[k])
    lp.set_max(lp[k],upper_bound)

  # objective function: in reality, we only want a feasible solution
  lp.set_objective(lp.sum([lp[k] for k in xrange(n)]))

  # add constraints relative to each choice of LT
  for i in xrange(len(G)):
    p = G[i].value()
    for m in p.monomials():
      if m != LTs[i]:
        vec_v1 = vector(LTs[i].exponents()[0])
        vec_v2 = vector(m.exponents()[0])
        constraint = tuple(vec_v1 - vec_v2)
        lp.add_constraint(lp.sum([constraint[k]*lp[k] for k in xrange(n)]),min=tolerance_cone)

  lp.solve()
  w = lp.get_values([lp[k] for k in xrange(n)])

  return w, new_vertices

@cython.profile(True)
cpdef tuple choose_cone_ordering(list G, list current_ordering, list constraints, MixedIntegerLinearProgram lp, set rejects, set bvs, int use_bvs, int use_dcs):
  r"""
  Choose an ordering in an unrestricted way but based on Perry's algorithm.
  Idea: reinserting previously computed polynomials, in order to choose new LMs for them.

  I think Perry's simplifications of constraints make so that this doesn't work
  evidence: too many cases where no constraints are added - i.e., only one LM was possible
  """

  #STEP 1: apply Perry's algorithm on the new polynomial
  cdef int i, j, k = len(G)
  cdef MPolynomialRing_libsingular R = G[0].value().parent()
  cdef int new_beginning = lp.number_of_constraints()
  cdef tuple result = choose_ordering_restricted(G, [ G[i].value().lm() for i in xrange(k-1)], k-1, current_ordering, lp, rejects, bvs, use_bvs, use_dcs, False)
  current_ordering = result[0]
  lp = result[1]
  bvs = result[2]
  cdef int new_end = lp.number_of_constraints()
  constraints.append((new_beginning, new_end))
  cdef MPolynomialRing_libsingular PR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(current_ordering))
  for i in xrange(len(G)):
    G[i].set_value(PR(G[i].value()))

  #STEP 2: build the lists of useful and useless LTs

  #We need to implement a criterion to decide which previous LTs are useless
  #TODO this is not very efficient. In practice, this could be computed once and updated on changes!
  cdef list useful_LTs = []
  cdef list useless_LTs = []
  for i in xrange(k):
    if any([ monomial_divides(G[j].value().lm(), G[i].value().lm()) for j in xrange(k) if j != i ]):
      useless_LTs.append(i) #lm of G[i] is useless
    else: #lm of G[i] is useful, add to list
      useful_LTs.append(G[i].value().lm())

  #print useless_LTs

  #STEP 3: Choose an old polynomial and reinsert it using Caboara and Perry's restricted criterion
  cdef int reinsert
  if len(useless_LTs) == 0:
    return result + (constraints, False)
  else:
    #reinsert = useless_LTs[randint(0, len(useless_LTs)-1)] #TODO This suffices as criterion for now, we can find a better one latter
    reinsert = randint(0, len(G)-2)

  #We need to keep a list of which constraints in lp correspond to which polynomials.
  #Then we remove relevant constraints from lp when a polynomial is removed and reinserted.

  #constraints is a list of tuples, each tuple with the beginning and end indices of constraints
  beginning, end = constraints[reinsert]
  print constraints
  if beginning != end:
    lp.remove_constraints(xrange(beginning, end))

    #Reindex stuff in constraints because some were removed
    for i in xrange(reinsert+1, len(constraints)):
      beg_i, end_i = constraints[i]
      constraints[i] = (beg_i - (end - beginning), end_i - (end - beginning))
  #Remove element from G and constraints list
  cdef clothed_polynomial g = G.pop(reinsert)
  constraints.pop(reinsert)

  print "showing g", g.value().lm()
  G.append(g)
  new_beginning = lp.number_of_constraints()
  result = choose_ordering_restricted(G, [G[i].value().lm() for i in xrange(k-1)], k-1, current_ordering, lp, rejects, bvs, use_bvs, use_dcs, False)
  new_end = lp.number_of_constraints()
  constraints.append((new_beginning, new_end))

  PR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(result[0]))

  return result + (constraints, PR(G[k-1].value()).lm() != G[k-1].value().lm())


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

# imports

from copy import copy
from random import randint, choice

import cython
from cython.parallel cimport prange

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix

from sage.libs.singular.decl cimport p_DivisibleBy

from sage.misc.misc_c import prod

#from sage.numerical.mip import Sum
from sage.numerical.mip cimport MixedIntegerLinearProgram

from sage.rings.infinity import Infinity

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular

from sage.rings.integer_ring import IntegerRing

from sage.modules.free_module_element import vector
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.matrix.matrix_real_double_dense cimport Matrix_real_double_dense
from sage.modules.vector_real_double_dense cimport Vector_real_double_dense
from sage.rings.real_double import RDF
from sage.geometry.polyhedron.constructor import Polyhedron

from sage.numerical.backends.glpk_backend cimport GLPKBackend

from sage.functions.other import floor, ceil

# globals, until/unless I make a class out of this

ZZ = IntegerRing()
# record keeping
cdef int rejections, monomials_eliminated, number_of_programs_created
# options
cdef int sugar_type
# tweaking the linear programs
cdef double tolerance_cone = 0.01
cdef double upper_bound = 100
cdef double upper_bound_delta = 100

# memory management courtesy C

cdef extern from "stdlib.h":
  ctypedef unsigned int size_t
  void *malloc(size_t size)
  void free(void *ptr)

# types

cdef class clothed_polynomial:

  r"""
  We use clothed polynomials to store information about polynomials,
  and to make it easier to change the ordering.
  """

  # internal data

  cdef MPolynomial_libsingular f # the polynomial that we have clothed
  cdef int sugar # the sugar of this polynomial, as computed by us

  def __init__(self, MPolynomial_libsingular f):
    r"""
    Initializes ``self.f`` to ``f`` and ``self.sugar`` to the computed sugar of ``f``.
    """
    global sugar_type
    cdef list exp = f.exponents(as_ETuples=False)
    cdef tuple tup
    cdef int d = 0

    self.f = f
    if sugar_type == 0:
      for tup in exp:
        d = max(d,sum(tup))
        self.sugar = d
    elif sugar_type == 1:
      self.sugar = f.degree()

  def __repr__(self): return self.f.__repr__()

  cdef is_equal(self, clothed_polynomial other): return self == other

  # methods related to the polynomial

  cpdef MPolynomial_libsingular value(self): return self.f

  cpdef set_value(self, MPolynomial_libsingular f): self.f = f

  cpdef lm(self): return self.f.lm()

  # methods related to sugar

  cpdef set_sugar(self, int s): self.sugar = s

  cpdef int get_sugar(self): return self.sugar

# utility functions

@cython.profile(True)
cpdef int monomial_divides(MPolynomial_libsingular t, MPolynomial_libsingular u):
  r"""
    I use this as a way to get around sage's ring method for divisibility,
    which includes a lot of error checking that, frankly, I don't need,
    and slows things down TREMENDOUSLY.
  """
  return p_DivisibleBy(t._poly, u._poly, t._parent_ring)

@cython.profile(True)
cpdef int indivisible(tuple t, tuple u):
  r"""
    Determines whether the tuple t represents a monomial indivisible
    by the monomial represented by the tuple u.
  """
  cdef int i = 0
  cdef int divisible = 1
  while divisible and i < len(t):
    if t[i] < u[i]: divisible = 0
    i += 1
  return not divisible

@cython.profile(True)
cpdef MixedIntegerLinearProgram new_linear_program(MixedIntegerLinearProgram lp = None):
  r"""
    This tracks the number of linear programs created, and initializes them
    with a common template.
  """
  global number_of_programs_created
  number_of_programs_created += 1
  if lp == None:
    return MixedIntegerLinearProgram(check_redundant=True, solver="GLPK", maximization=False)
  else: return copy(lp)

@cython.profile(True)
cpdef tuple monitor_lts(list G, list LTs, list new_ordering):
  r"""
    Checks whether new_ordering changes the leading terms of G.

    INPUT:

    - ``G`` -- current basis
    - ``LTs`` -- leading terms of elements of G according to previous ordering
    - ``new_ordering`` -- new ordering

    OUTPUT:

    - 1 if polynomials whose leading monomials change; 0 otherwise
    - list of lists; if first result is nonzero, then each list in this list
      contains the exponent vectors of monomials that weigh the same or more
      than the correct leading monomial
  """
  cdef int i, j, k, n, changes
  cdef double tw, uw # monomial weights
  cdef tuple texp, uexp # monomial exponents
  cdef MPolynomial_libsingular g
  cdef list U, result, current
    # U is the set of monomials of the polynomial g currently under examination
    # current is the list of monomials of g that weigh more than the old polynomial
    # result is collection of all current's

  # setup

  n = len(new_ordering)
  result = list()
  changes = 0

  # main loop

  # check each leading monomial
  for i in xrange(len(LTs)):

    texp = LTs[i].exponents(as_ETuples=False)[0]
    tw = sum([texp[k]*new_ordering[k] for k in xrange(n)])
    current = list()
    result.append(current)
    g = G[i].value(); U = g.monomials()

    # heaviest entry in U should be t
    for u in U:

      uexp = u.exponents(as_ETuples=False)[0]

      if uexp != texp:

        uw = sum([uexp[k]*new_ordering[k] for k in xrange(n)])

        if uw >= tw:

          changes = 1
          current.append(uexp)

  #print result
  return changes, result

@cython.profile(True)
cpdef set boundary_vectors(MixedIntegerLinearProgram lp, int n):
  r"""
    Finds a boundary vector for an admissible ordering defined by ``lp``,
    so that a vector satisfying ``lp`` passes between these vectors,
    WITH HIGH PROBABILITY. By this, I mean that we may lose some vectors,
    but as long as ``lp`` is defined well, this should be rare in general.
    It is ESSENTIAL that the program have at least one solution, because
    we WILL NOT check that here.

    INPUT:

    - ``lp`` -- a linear program corresponding to an admissible ordering
    - ``n`` -- number of variables in the program

    OUTPUT:

    A tuple containing a list of corner vectors that approximate the region

    ALGORITHM:

      #. find the minimum feasible degree d;
         that is, there exists a feasible solution (x1,...,xn) such that
         x1 + ... + xn = d
      #. create a cross-section of the solution cone
         by intersecting with the hyperplane x1 + ... + xn = d + 1
      #. for each variable x, compute the vectors that maximize and minimize x
         in this cross-section
  """
  cdef int i, j, k # counters
  cdef int start_constraints, end_constraints # used for deleting bad constraints
  cdef float sol_degree # degree of solution
  cdef tuple result, solution_vector
  cdef set boundaries
  cdef list lp_sol

  #print "in boundary vectors"
  np = lp
  start_constraints = np.number_of_constraints()

  # first we compute the minimum feasible degree
  # this first step should add zero overhead
  np.solve()
  lp_sol = np.get_values([np[k] for k in xrange(n)])
  sol_degree = sum(lp_sol)

  # now create the cross-section
  np.add_constraint(np.sum([np[k] for k in xrange(n)]), min=sol_degree + 1)
  np.add_constraint(np.sum([np[k] for k in xrange(n)]), max=sol_degree + 1)
  end_constraints = start_constraints + 2

  # find corners where variables are maximized and minimized
  boundaries = set()

  for k in xrange(n):

    # lp is a minimization problem, so this minimizes xk
    np.set_objective(np[k])
    #print "variable", k, "min"
    np.solve()
    boundaries.add(tuple(np.get_values([np[k] for k in xrange(n)])))
    # lp is a minimization problem, so this maximizes xk
    np.set_objective(-np[k])
    #print "variable", k, "max"
    np.solve()
    boundaries.add(tuple(np.get_values([np[k] for k in xrange(n)])))

  # now remove the cross-section
  np.remove_constraints((start_constraints,start_constraints+1))
  np.set_objective(np.sum([np[k] for k in xrange(n)]))
  np.solve()

  print "boundaries", boundaries
  #print "leaving boundary vectors"
  return boundaries

@cython.profile(True)
cpdef int solve_real(MixedIntegerLinearProgram lp, int n):
  r"""
  Set the linear program ``lp`` to treat each of the ``n`` variables
  as a real variable, not an integer variable. Then, solve the program.
  """
  cdef int k
  global failed_systems

  for k in xrange(n): lp.set_real(lp[k])

  try: lp.solve()
  except:
    failed_systems += 1
    return 0

  return 1

@cython.profile(True)
cpdef int solve_integer(MixedIntegerLinearProgram lp, int n):
  r"""
  Set the linear program ``lp`` to treat each of the ``n`` variables
  as an integer variable, not a real variable. Then, solve the program.

  A difficulty in using glpk is that it requires an upper bound in order to find
  an integer solution in a decent amount of time.
  As a result, we have set up our linear program with such an upper bound.
  However, while the real solution can be found with even a fairly small upper bound,
  the integer solution will at times exceed the previously set upper bound.
  This poses a challenge.

  This function is not called if there is not a real solution,
  and the structure of our feasible regions implies that
  if there is a real solution, then there is also an integer solution.
  If no solution is found here, then, the problem must be that our upper bound
  is too low.
  """
  global upper_bound, upper_bound_delta, failed_systems
  cdef int k, passed
  cdef double old_upper_bound

  for k in xrange(n):
    lp.set_integer(lp[k])
    lp.set_max(lp[k], upper_bound)

  # try to get solution
  passed = False
  old_upper_bound = upper_bound

  while not passed:

    # first find ordering
    print "obtaining ordering below", upper_bound

    try:

      lp.solve()
      passed = True

    except:

      # upper bound too low, so raise it
      print "failed"
      upper_bound += upper_bound_delta

      if upper_bound > 128000: # even I give up after a certain point -- need a better solver?

        print "returning no solution", upper_bound
        for k in xrange(n): lp.set_real(lp[k])
        upper_bound = old_upper_bound
        failed_systems += 1
        return 0

      for k in xrange(n): lp.set_max(lp[k], upper_bound)

  return 1

@cython.profile(True)
cpdef tuple feasible(int i, list CMs, MixedIntegerLinearProgram olp, set rejects, list G, list LTs, int use_rejects):
  r"""
    Determines if the ``i``th monomial in ``CMs`` is a feasible leading monomial
    of the polynomial whose compatible monomials are listed in CMs,
    consistent with the linear program olp.

    INPUT:

    - ``i`` -- an integer (0 <= i < len(M))
    - ``CMs`` -- monomials judged to be compatible with previous choices
    - ``olp`` -- a linear program for solutions
    - ``rejects`` -- a set of linear constraints that are known to be incompatible
      with ``olp``
    - ``G`` -- the current basis of the ideal
    - ``LTs`` -- the leading terms of ``G``
    - `use_rejects` -- whether to use rejected programs (disjoint cones) to avoid useless systems

    OUTPUT:

    A tuple containing nothing, if there is no solution. Otherwise, it contains:

      * ``i``,
      * a vector that solves ``lp``, and
      * ``lp``, a linear program that extends ``olp`` and whose solution vector
        v satisfies v.CMs[i] > v.CMs[j] for all j=/=i.
        (Here, a.b indicates the dot product.)

    ALGORITHM:

    It constructs two linear systems, both of the form

      #. `x_k \geq \epsilon` for each `k`
      #. `\sum((alpha_k - beta_k)*x_k) \geq \epsilon` where

        #. `alpha_k` is the degree of `t` in `x_k` and
           `beta_k` is the degree of `u` from `CMs[i]` in `x_k` and
        #. epsilon is a constant used to minimize floating point error;
           its value is determined by the global `tolerance_cone`

    The first system merely checks against elements of CMs; that is,
    against the current polynomial.
    The second combines the first system with ``olp``; that is,
    against all the previous polynomials, as well.

  """
  # the following line is needed to control the solving process
  import sage.numerical.backends.glpk_backend as glpk_backend

  global rejections, monomials_eliminated, tolerance_cone, upper_bound, failed_systems

  cdef int j, k, l # counters
  cdef int passed, lms_changed # signals
  cdef int n, start_constraints, end_constraints, number_of_constraints # measures
  cdef float sol_degree, a
  cdef tuple t, u, v, old_lm, new_lm # exponent vectors
  cdef tuple constraint
  cdef list changed_lms, changed_lm # used for when we change the ordering, rather than refine it

  cdef set new_constraints = set()
  cdef set new_rejects = set()

  t = <tuple>CMs[i]
  #print "testing", t
  n = len(t)
  cdef tuple result = tuple()

  # STEP 1: check if solution exists merely for this polynomial

  # set up solver to solve LP relaxation first
  # (GLPK chokes on no integer solution)
  cdef MixedIntegerLinearProgram lp = new_linear_program()
  lp.solver_parameter(glpk_backend.glp_simplex_or_intopt, glpk_backend.glp_simplex_then_intopt)

  # xi >= epsilon
  for k in xrange(n): lp.add_constraint(lp[k],min=tolerance_cone)
  # minimize x1 + ... + xn
  lp.set_objective(lp.sum([lp[k] for k in xrange(n)]))
  #print "initial constraints set"

  # add constraints
  cdef int m = len(CMs)
  #print m, "new constraints", CMs

  #loop through all exponent vectors of compatible monomials
  for j in xrange(m):

    u = <tuple>CMs[j] # get current exponent vector

    if t != u: # don't check against oneself!

      constraint = tuple([float(t[k]-u[k]) for k in xrange(n)])
      new_constraints.add(constraint)

      # check whether adding constraint to this system is known to be incompatible
      if use_rejects:
        #print "checking rejects", 
        for c in (rejects):
  
          if c.issubset(new_constraints):
  
            #print "rejected", new_constraints, c
            #print "rejected!!!"
            rejections += 1
            return result

      #print "checked rejects"
      lp.add_constraint(lp.sum([constraint[k]*lp[k] for k in xrange(n)]),min=tolerance_cone)

      try:

        #print "LP relaxation"
        lp.solve()

      except: # unable to solve; monomial cannot lead even this polynomial

        failed_systems += 1
        rejects.add(frozenset(new_constraints)) # remember this failed system
        return result

  # we have a solution to the LP relaxation
  #print "have solution"

  # STEP 2: check if system is consistent with previous systems

  # first, solve the LP relaxation
  #print "copying program"
  np = olp; start_constraints = end_constraints = np.number_of_constraints()
  #print "copied old program, which had", start_constraints, "constraints"

  # add new constraints to old constraints
  for constraint in new_constraints:

    #print t, constraint
    number_of_constraints = np.number_of_constraints()
    new_rejects.add(constraint)
    np.add_constraint(np.sum([constraint[k]*np[k] for k in xrange(n)]), min=tolerance_cone)
    # check of size is necessary because we prevent the addition of redundant constraints
    if np.number_of_constraints() > number_of_constraints: end_constraints += 1

  #print "now have", np.number_of_constraints(), "constraints"

  # moving this inside the loop might make rejection-checking more efficient,
  # but slows the process greatly for large systems
  try:

    # a lot of debugging information for when I mess things up...
    #print "trying to solve with old constraints"
    #np.show()
    np.solve()
    #print "solved"

  except: # failed to find a solution; monomial is incompatible w/previous

    #print "not solved"
    rejects.add(frozenset(new_rejects)) # remember this failed system
    #print "killing new constraints in program, which has", np.number_of_constraints(), "constraints"
    np.remove_constraints(range(start_constraints,end_constraints))
    failed_systems += 1
    return result

  # now, try integer optimization
  #print "have real solution", np.get_values([np[k] for k in xrange(n)])
  #print "copied constraints"
  cdef list sol

  # STEP 3: obtain weight vector, which must have integer solutions

  upper_bound = round(sum(np.get_values([np[k] for k in xrange(n)]))) + 1

  if not solve_integer(np, n): return False

  # make sure older LTs have not changed
  sol = np.get_values([np[k] for k in xrange(n)])
  lms_changed = 1

  while lms_changed != 0:

    #print "checking lms"
    lms_changed, changed_lms = monitor_lts(G, LTs, sol)
    resolve = False

    for j in xrange(len(changed_lms)):

      if len(changed_lms[j]) != 0: # some monomial changed :-(

        #print "adding constraints for", j
        resolve = True
        u = changed_lms[j][0]
        #print u
        v = LTs[j].exponents(as_ETuples=False)[0]
        constraint = tuple([float(v[k] - u[k]) for k in xrange(n)])
        new_constraints.add(constraint)
        number_of_constraints = np.number_of_constraints()
        np.add_constraint(np.sum([constraint[k]*np[k] for k in xrange(n)]), min=tolerance_cone)
        if np.number_of_constraints() > number_of_constraints: end_constraints += 1

    # if a monomial changed, we have to solve anew
    #print "resolve?", resolve
    if resolve:

      #print "resolving"
      if solve_real(np, n) and solve_integer(np, n):

        #print "found solutions"
        sol = np.get_values([np[k] for k in xrange(n)])

      else: # a-ha! this polynomial was actually incompatible!

        #print "failed to solve"
        np.remove_constraints(range(start_constraints, end_constraints))
        rejects.add(frozenset(new_constraints))
        return result

  #print "got it", sol
  # return to LP relaxation for future work
  for k in xrange(n):

    np.set_real(np[k])
    np.set_max(np[k],64000)

  # set up new solution
  result = (t,sol,np)
  print np.number_of_constraints(), "constraints"
  #print "returning solution"
  return result

@cython.profile(True)
cpdef list possible_lts(MPolynomial_libsingular f, set boundary_vectors, int use_boundary_vectors):
  r"""
    Identifies terms of ``f`` that could serve as a leading term,
    and are compatible with the boundary approximated by ``boundary_vectors``.

    INPUT:

    - ``f`` -- a polynomial
    - ``boundary_vectors`` -- a tuple of maximum and minimum values of the variables
      on the current feasible region, allowing us to approximate the corners of
      the solution cone
    - `use_boundary_vectors` -- whether to use the boundary vectors; setting this to
      `False` gives behavior similar to Caboara's original implementation

    OUTPUT:

    A list of tuples, representing the exponent vectors of the potential leading monomials.

    ALGORITHM:

      #. Let ``t`` be current leading monomial of ``f``.
      #. Let ``U`` be set of exponent vectors of monomials of ``f``, except ``t``.
      #. Let ``M`` be subset of ``U`` s.t. `u\in M` iff `c.u > c.t`
        for some `c\in boundary_vectors`. (here, c.u denote multiplication)
      #. Return ``M``.

  """
  #print "in possible_lts"
  cdef list M
  cdef MPolynomial_libsingular ux, vx # terms
  cdef int i, j # counters
  cdef int passes # signals
  cdef tuple ordering, u, v # ordering and exponent vectors
  global monomials_eliminated

  # identify current leading monomial, other monomials; obtain their exponent vectors

  cdef MPolynomialRing_libsingular R = f.parent()
  cdef int n = len(R.gens())
  cdef tuple t = f.lm().exponents(as_ETuples=False)[0]
  cdef list U = f.monomials()
  U.remove(f.lm())

  # determine which elements of U might outweigh t somewhere within the current solution cone
  # if there is no solution cone yet (first polynomial), pass them all

  # no point in checking anything if we have no boundaries yet...
  # todo: initialize boundary_vectors to obvious boundaries (e_i for i=1,...,n) on first polynomial
  # so that we can remove this check
  if use_boundary_vectors and boundary_vectors != None:

    M = list()

    for j in xrange(len(U)):

      ux = U[j]; u = ux.exponents(as_ETuples=False)[0]; passes = False

      # check whether u weighs more than t according to some boundary vector
      for ordering in boundary_vectors:

        # compare dot products
        if sum([u[k]*ordering[k] for k in xrange(n)]) > sum([t[k]*ordering[k] for k in xrange(n)]):

          M.append((u,ux))
          passes = True
          break

      if not passes: monomials_eliminated += 1

  else: M = [(ux.exponents(as_ETuples=False)[0], ux) for ux in U] # if no boundary_vectors, allow all monomials

  # don't forget to append t
  M.append((t,f.lm()))

  # remove monomials u divisible by some other monomial v -- this will happen rarely
  # when using boundary vectors, but could happen even then
  cdef list V = list()

  for (u,ux) in M:

    passes = True; i = 0

    while passes and i < len(M):

      v, vx = M[i]
      if vx != ux and monomial_divides(ux, vx): passes = False
      i += 1

    if passes: V.append((u,ux))

  print V
  #print len(M), "possible leading monomials"
  return V

cpdef int hs_heuristic(f, g):
  r"""
    Implements the Hilbert heuristic recommended by Caboara in his 1993 paper.
    It first compares the degrees of the Hilbert polynomials;
    if one polynomial's degree is smaller than the other,
    then it returns the difference b/w the degree for f's hp and g's.
    Otherwise, it compares the Hilbert series of the polynomials;
    if the two are equal, it returns 0; otherwise,
    it returns the trailing coefficient of the numerator of the difference b/w f's hs and g's.
    See Caboara's paper for details
  """

  if f[0].degree() == g[0].degree():

    if f[1] == g[1]: return 0
    else:
      C = (f[1]-g[1]).numerator().coefficients()
      return C[len(C)-1]

  else: return f[0].degree() - g[0].degree()

@cython.profile(True)
cpdef list sort_CLTs_by_Hilbert_heuristic(MPolynomialRing_libsingular R, list current_Ts, list CLTs):
  r"""
    Sorts the compatible leading monomials using the Hilbert heuristic
    recommended by Caboara in his 1993 paper.
    Preferable leading monomials appear earlier in the list.

    INPUTS:
    - `R` -- current ring (updated according to latest ordering!)
    - `current_Ts` -- current list of leading terms
    - `CLTs` -- compatible leading terms from which we want to select leading term
      of new polynomial
  """
  cdef tuple tup

  # create a list of tuples
  # for each leading term tup, we consider the monomial ideal created if tup
  # is chosen as the leading term.
  # the first entry is the tentative Hilbert polynomial
  # the second is the tentative Hilbert series
  # the third is tup itself (the compatible leading term)
  CLTs = [(R.ideal(current_Ts + [tup[1]]).hilbert_polynomial(),
             R.ideal(current_Ts + [tup[1]]).hilbert_series(),
             tup)
           for tup in CLTs]
  CLTs.sort(cmp=hs_heuristic) # sort according to hilbert heuristic
  #print CLTs

  return CLTs

cpdef list min_CLT_by_Hilbert_heuristic(MPolynomialRing_libsingular R, list CLTs):
  r"""
    Sorts the compatible leading monomials using the Hilbert heuristic
    recommended by Caboara in his 1993 paper.
    Preferable leading monomials appear earlier in the list.

    INPUTS:
    - `R` -- current ring (updated according to latest ordering!)
    - `CLTs` -- compatible leading terms from which we want to select leading term
      of new polynomial

    OUTPUTS:
    - the preferred list of leading monomials from CLTs
  """

  cdef list leads #lists of leading monomials in CLTs

  CLTs = [(R.ideal(leads).hilbert_polynomial(),
           R.ideal(leads).hilbert_series(),
           leads)
          for leads in CLTs]
  CLTs.sort(cmp=hs_heuristic)
  #TODO this could be optimized, we are sorting just to find the minimum...
  return CLTs[0][2]

cpdef list min_weights_by_Hilbert_heuristic(MPolynomialRing_libsingular R, list CLTs):
  r"""
    Sorts the compatible leading monomials using the Hilbert heuristic
    recommended by Caboara in his 1993 paper.
    Preferable leading monomials appear earlier in the list.

    INPUTS:
    - `R` -- current ring (updated according to latest ordering!)
    - `current_Ts` -- current list of leading terms
    - `CLTs` -- compatible leading terms from which we want to select leading term
      of new polynomial

    OUTPUTS:
    - the preferred weight vector from CLTs
  """

  cdef tuple leads #lists of leading monomials in CLTs

  CLTs = [(R.ideal(leads[1]).hilbert_polynomial(),
           R.ideal(leads[1]).hilbert_series(),
           leads[0])
          for leads in CLTs]
  CLTs.sort(cmp=hs_heuristic)
  #TODO this could be optimized, we are sorting just to find the minimum...
  return CLTs[0][2]

cpdef GLPKBackend make_solver():
  r"""
  Creates an empty model in a GLPK backend solver.

  OUTPUTS:
  - a GLPK backend
  """

  from sage.numerical.backends.generic_backend import get_solver
  import sage.numerical.backends.glpk_backend as backend
  lp = get_solver(solver="GLPK")
  lp.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)

  return lp

cpdef void append_linear_program(GLPKBackend glpk, MPolynomial_libsingular p):
  r"""
  Appends constraints and variables of the Newton polyhedron of `p` to the current
  linear program in `glpk`.

  INPUTS:

  - `glpk` -- current representation of the linear programming model in GLPK
  - `p` -- polynomial whose affine Newton Polyhedron will be added to the linear programming model
  """

  cdef int n = p.parent().ngens()
  cdef tuple a
  cdef list variables, coefs

  NP = p.newton_polytope() + Polyhedron(rays=(-identity_matrix(n)).rows())
  lp = NP.to_linear_program(solver="GLPK")
  n = glpk.ncols()
  glpk.add_variables(lp.number_of_variables())

  for lb, a, ub in lp.constraints():
    variables, coefs = a
    variables = [ i + n for i in variables ]
    glpk.add_linear_constraint(list(zip(variables, coefs)), lb, ub)

  return

@cython.profile(True)
cpdef list weight_vector(GLPKBackend lp, int n):
  r"""
  Returns the weight vector currently used as objective function in the linear
  programming model lp.

  INPUTS:

  - `lp` -- a GLPK linear programming model

  OUTPUTS:

  - a list representing a weight vector
  """
  cdef list w = []
  cdef int i
  for i in xrange(n):
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
cpdef list sensitivity(GLPKBackend lp, int n, int k):
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

  assert all([ z >= 0 for z in zN ]), "Current solution is not optimal?"

  cdef list DcB_l = []
  for i in basic:
    if i >= m and (i - m) % n == change:
      DcB_l.append(i)

  cdef Vector_real_double_dense DzN = sum([ tableau_row(lp, i, nonbasic) for i in DcB_l ]) - DcN

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

  print lower, upper
  assert upper >= 0 and lower <= 0, "Sensibility range is inconsistent"

  #Change the objective function according to range found
  cdef float increment

  if lower == float("-inf") and upper == float("inf"):
    return weight_vector(lp, n)
  if upper == float("inf"):
    increment = ceil(lower) - 1
    if lp.objective_coefficient(change) + increment <= 0:
      return weight_vector(lp, n)
  elif lower == float("-inf"):
    increment = floor(upper) + 1
  else:
    #Lower coefficient, when possible, 50% of the time
    increment = ceil(lower) - 1
    if lp.objective_coefficient(change) + increment <= 0:
      increment = floor(upper) + 1
    elif randint(0, 1):
      increment = floor(upper) + 1

  cdef float old_value = lp.objective_coefficient(change)
  for i in xrange(k):
    lp.objective_coefficient(change + i * n, old_value + increment)
  return weight_vector(lp, n)

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

#TODO why does the number of iterations affect performance so much?
@cython.profile(True)
cpdef list choose_simplex_ordering(list G, list current_ordering, GLPKBackend lp, int iterations = 5):
  r"""

  INPUTS:

  - `G` -- the current system of generators
  - `current_ordering` -- the current weight ordering
  - `lp` -- previous GLPK linear programming model
  - `iterations` -- number of neighbors to visit

  OUTPUTS:

  - a list of weights representing a monomial order
  """
  cdef MPolynomialRing_libsingular R = G[0].value().parent()
  cdef MPolynomialRing_libsingular newR
  cdef int k = len(G)
  cdef int n = R.ngens()
  cdef int i, j, it = 0
  cdef list CLTs, LTs, oldLTs, w, best_w

  #Initial random ordering
  w = [ randint(1, 10000) for i in xrange(n) ]
  best_w = w

  #Transform last element of G to linear program, set objective function given by w and solve
  append_linear_program(lp, G[len(G)-1].value())
  lp.set_objective(w * k)
  lp.solve()

  #Get current LTs to compare with Hilbert heuristic
  newR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(w))
  LTs = find_monomials(lp, newR, k)
  oldLTs = LTs
  CLTs = [ (newR.ideal(LTs).hilbert_polynomial(), newR.ideal(LTs).hilbert_series(), w ) ]

  #Do sensitivity analysis to get neighbor, compare
  while it < iterations:
    #print "before"
    #L = []
    #for i in xrange(k):
    #  L2 = []
    #  for j in xrange(n):
    #    L2.append(lp.get_variable_value(j + i * n))
    #  L.append(L2)
    #print L
    w = sensitivity(lp, n, k)
    lp.solve()
    #print "after"
    #L = []
    #for i in xrange(k):
    #  L2 = []
    #  for j in xrange(n):
    #    L2.append(lp.get_variable_value(j + i * n))
    #  L.append(L2)
    #print L
    newR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(w))
    LTs = find_monomials(lp, newR, k)
    #print w, best_w
    print [LTs[i] == oldLTs[i] for i in xrange(len(LTs))].count(False), len(LTs)
    CLTs.append((newR.ideal(LTs).hilbert_polynomial(), newR.ideal(LTs).hilbert_series(), w))
    CLTs.sort(cmp=hs_heuristic)
    best_w = CLTs[0][2] #Take first improvement
    if best_w == w:
      oldLTs = LTs
    CLTs = CLTs[:1]
    it += 1
    lp.set_objective(best_w * k)
    lp.solve()

  #Compare with current_ordering - keep the current one if they tie!
  newR = PolynomialRing(R.base_ring(), R.gens(), order=create_order(current_ordering))
  LTs = [ newR(G[k].value()).lm() for k in xrange(len(G)) ]
  CLTs.insert(0, (newR.ideal(LTs).hilbert_polynomial(), newR.ideal(LTs).hilbert_series(), current_ordering))
  CLTs.sort(cmp=hs_heuristic)
  best_w = CLTs[0][2]
  lp.set_objective(best_w * k)
  lp.solve()

  return best_w

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
cpdef tuple choose_ordering_restricted(list G, list current_Ts, int mold, list current_ordering, MixedIntegerLinearProgram lp, set rejects, set bvs, int use_bvs, int use_dcs):
  r"""
    Chooses a weight vector for a term ordering for the basis ``G`` that refines
    the weight vector that solves ``lp``.

    INPUTS:

    - ``G`` --  a basis of a polynomial ideal
    - ``current_Ts`` -- leading monomials of *some* elements of ``G``,
      according to the current ordering
    - ``mold`` -- the number of elements of ``current_Ts``
    - ``lp`` -- a linear program whose solutions dictate the choice ``current_Ts``
    - ``rejects`` -- record of linear programs rejected previously
    - ``bvs`` -- approximation to feasible region, in the form of
      maximum and minimum values of a cross-section of the feasible region
    - `use_bvs` -- whether to use boundary vectors; setting this and `use_dcs` to False
      gives us behavior similar to Caboara's original implementation
    - `use_dcs` -- whether to use disjoint cones; setting this and `use_bvs` to False
      gives us behavior similar to Caboara's original implementation

    OUTPUTS:

    - a weighted ordering that solves...
    - a linear program that extends ``lp``, approximated by...
    - a list of vectors

    ALGORITHM:

    #. Of the possible leading monomials of each new `g\in G` that are compatible with ``lp``,
    #. determine which combinations are consistent with each other, then
    #. identify one which we think is a good choice.
  """
  #print "in optimal ordering"
  cdef int i, j, k # counters

  cdef MPolynomial_libsingular t, g # term, polynomial
  cdef MPolynomialRing_libsingular R = G[0].value().parent() # current ring

  # measures
  cdef int m = len(G)
  cdef int n = len(R.gens())

  #signals
  cdef int found, passes

  # contain leading terms, exponent vectors, weight ordering
  cdef list CLTs, LTs, LTups, w
  # elements of previous
  cdef tuple can_work, tup, changed_lms

  #print mold, m
  w = current_ordering

  for i in xrange(mold,m): # treat each new polynomial

    #print "finding lts for", i
    g = G[i].value()
    #print len(g.monomials()), "monomials to start with"
    CLTs = possible_lts(g, bvs, use_bvs)
    #print len(CLTs), "compatible leading monomials"
    #print CLTs

    #print("before sort")
    # use Caboara's Hilbert heuristic
    CLTs = sort_CLTs_by_Hilbert_heuristic(R, current_Ts, CLTs)
    #print CLTs
    # discard unnecessary information -- need to see why I have this information in the first place
    CLTs = [tup[len(tup)-1] for tup in CLTs]
    # extract the leading tuples
    LTups = [tup[0] for tup in CLTs]

    # find a compatible leading monomial that actually works
    # use fact that preferred monomials appear earlier in list
    found = False
    j = 0

    while not found and CLTs[j][1] != g.lm():

      print "testing", CLTs[j], "against", g.lm()
      #print j, LTups, CLTs
      can_work = feasible(j, LTups, lp, rejects, G, current_Ts, use_dcs)
      #print can_work

      if len(can_work) > 0: # this means we found a solution

        w = can_work[1]
        print CLTs[j][1], "worked! with", can_work
        found = True

      else: # monomial turned out to be incompatible after all

        print CLTs[j][1], "did not work!"
        j += 1
    # now that we've found one, use it
    #print "success with monomial", j, CLTs[j], LTups[j]
    # no point in finding boundary vectors if g has the same leading term as before
    if CLTs[j][1] == g.lm(): return current_ordering, lp, bvs
    t = <MPolynomial_libsingular>CLTs[j][1]
    current_Ts.append(t) # hmm
    G[i].set_value(G[i].value() * G[i].value().coefficient(t)**(-1)) # make polynomial monic
    lp = can_work[2] # get new linear program
    #lp.show()
    bvs = boundary_vectors(lp, n) # compute new boundary vectors
    w = can_work[1] # get new weight vector

  #print "We have", len(rejects), "rejects"
  return w, lp, bvs

@cython.profile(True)
cpdef clothed_polynomial spoly(tuple Pd, list generators):
  r"""
    Computes an s-polynomial indexed by ``Pd``.

    INPUT:

    - ``Pd`` -- tuple (f,g); we wish to compute spoly(`f,g`) where f and g are clothed
    - `generators` -- the generators of the ideal
  """
  global sugar_type
  # counters
  cdef int k
  # measures
  cdef int new_sugar
  # polynomials, clothed and otherwise
  cdef MPolynomialRing_libsingular R
  cdef MPolynomial_libsingular f, g, tf, tg, tfg
  cdef clothed_polynomial cf, cg, cs

  # get polynomials
  cf = Pd[0]; cg = Pd[1]
  f = cf.value(); g = cg.value()

  if g != 0: # is f a generator?

    # no; get leading monomials of both f and g, then lcm
    tf = f.lm(); tg = g.lm()
    tfg = tf.lcm(tg)
    #print "building", (tf, tg, tfg)
    R = tfg.parent()
    s = R.monomial_quotient(tfg,tf)*f - R.monomial_quotient(tfg,tg)*g

    if sugar_type == 0: # standard sugar: add exponents
      new_sugar = max(cf.get_sugar() + sum(R.monomial_quotient(tfg,tf).exponents(as_ETuples=False)[0]),
        cg.get_sugar() + sum(R.monomial_quotient(tfg,tg).exponents(as_ETuples=False)[0]))

    elif sugar_type == 1: # weighted sugar: add degrees (based on ordering)
      new_sugar = max(cf.get_sugar() + R.monomial_quotient(tfg,tf).degree(),
        cg.get_sugar() + R.monomial_quotient(tfg,tg).degree())

  else: # f is a generator

    # s-polynomial IS f
    s = f; new_sugar = cf.get_sugar()
    k = 0
    print "building", f.lm()

    while k < len(generators):
      if generators[k] == cf: generators.pop(k)
      k += 1

  #print s
  cs = clothed_polynomial(s); cs.set_sugar(new_sugar)
  return cs

@cython.profile(True)
cpdef clothed_polynomial reduce_polynomial_with_sugar(clothed_polynomial s, list G):
  r"""
    Value of ``cf`` changes to reduced polynomial, and new sugar is computed.
    based on algorithm in Cox, Little, and O'Shea, and paper on Sugar by Giovini et al.

    INPUTS:
    - `s` -- s-polynomial to reduce
    - `G` -- list of reducers
  """
  global sugar_type
  # counters
  cdef int i, m
  # measures
  cdef int d
  # signals
  cdef int divided
  # polynomials, clothed and otherwise
  cdef clothed_polynomial cg
  cdef MPolynomial_libsingular f, g, r
  cdef MPolynomialRing_libsingular R
  cdef MPolynomial_libsingular t

  if s.value() == 0 or len(G) == 0: return s # no sense wasting time
  f = s.value()
  #TO1 = f.parent().term_order()
  d = s.get_sugar()
  R = f.parent()
  m = len(G)
  r = R(0)

  while f != 0: # continue reducing until nothing is left

    i = 0; divided = False

    # print "reducing with", f.lm()
    while f != 0 and i < m and not divided: # look for a monomial to reduce

      g = G[i].value()

      if monomial_divides(g.lm(), f.lm()):

        t = R.monomial_quotient(f.lm(),g.lm())
        f -= f.lc()/g.lc() * t * g

        if sugar_type == 0: # standard sugar: use exponents
          d = max(d, G[i].get_sugar() + sum(t.exponents(as_ETuples=False)[0]))

        elif sugar_type == 1: # weighted sugar: use degree (determined by ordering)
          d = max(d, G[i].get_sugar() + t.degree())

        divided = True
        i = 0

      else: i += 1

    if not divided: # did not divide; add to remainder
      r += f.lc() * f.lm()
      f -= f.lc() * f.lm()

  #print r
  # completed reduction; clean up
  if r != 0: r *= r.lc()**(-1) # make monic
  s.set_value(r); s.set_sugar(d)
  return s

@cython.profile(True)
cpdef clothed_polynomial reduce_poly(clothed_polynomial s, list G):
  r"""
    Reduces the s-polynomial `s`` modulo the polynomials ``G``
  """
  # polynomials
  cdef MPolynomial_libsingular r

  if s.value() != 0: # no point wasting time
    r = s.value().reduce(G)
    s.set_value(r)

  return s

@cython.profile(True)
cpdef int sug_of_critical_pair(tuple pair):
  r"""
    Compute the sugar of a critical pair.
  """
  global sugar_type
  # measures
  cdef int sf, sg, sfg, su, sv
  # polynomials, clothed and otherwise
  cdef clothed_polynomial cf, cg
  cdef MPolynomial_libsingular f, g, tf, tg, tfg, u, v
  # base ring
  cdef MPolynomialRing_libsingular R

  # initialization
  cf, cg = pair
  f, g = cf.value(), cg.value()
  tf = f.lm(); tg = g.lm()
  R = g.parent()

  # watch for a generator
  if f == 0:
    sfg = cg.get_sugar(); tfg = tg

  elif g == 0:
    sfg = cf.get_sugar(); tfg = tf

  else: # compute from how s-polynomial is constructed

    tfg = tf.lcm(tg)
    u = R.monomial_quotient(tfg, tf); v = R.monomial_quotient(tfg, tg)
    sf = cf.get_sugar(); sg = cg.get_sugar()

    if sugar_type == 0: # standard sugar: use exponents
      su = sum(u.exponents(as_ETuples=False)[0]); sv = sum(v.exponents(as_ETuples=False)[0])

    elif sugar_type == 1: # weighted sugar: use degree, determined by ordering
      su = u.degree(); sv = v.degree()

    sfg = max(sf + su, sg + sv)

  return sfg

@cython.profile(True)
cpdef lcm_of_critical_pair(tuple pair):
  r"""
    Compute the lcm of a critical pair.
  """
  # polynomials, clothed and otherwise
  cdef clothed_polynomial cf, cg
  cdef MPolynomial_libsingular f, g, tf, tg, tfg

  # initialization
  cf, cg = pair[0], pair[1]
  f, g = cf.value(), cg.value()
  tf = f.lm(); tg = g.lm()
  R = g.parent()

  if f == 0: tfg = tg
  elif g == 0: tfg = tf
  else: tfg = tf.lcm(tg)

  return tfg

@cython.profile(True)
cpdef deg_of_critical_pair(tuple pair):
  r"""
    Compute the exponent degree of a critical pair, based on the lcm.
  """
  # measures
  cdef int sfg
  # polynomials, clothed and otherwise
  cdef clothed_polynomial cf, cg
  cdef MPolynomial_libsingular f, g, tf, tg, tfg

  # initialization
  cf, cg = pair[0], pair[1]
  f, g = cf.value(), cg.value()
  tf = f.lm(); tg = g.lm()

  if f == 0: tfg = tg
  elif g == 0: tfg = tf
  else: tfg = tf.lcm(tg)

  sfg = sum(tfg.exponents(as_ETuples=False)[0])
  return sfg

# the next three functions are used for sorting the critical pairs

@cython.profile(True)
cpdef last_element(tuple p): return p[len(p)-1] # b/c cython doesn't allow lambda experessions, I think

@cython.profile(True)
cpdef lcm_then_last_element(tuple p): return (lcm_of_critical_pair(p), p[len(p)-1])

@cython.profile(True)
cpdef last_element_then_lcm(tuple p): return (p[len(p)-1], lcm_of_critical_pair(p))

@cython.profile(True)
cpdef gm_update(MPolynomialRing_libsingular R, list P, list G, list T, strategy):
  r"""
  The Gebauer-Moeller algorithm.
  
  INPUTS:
  
  - ``R`` -- the current ring (used for efficient division)
  - ``P`` -- list of critical pairs
  - ``G`` -- current basis (to discard redundant polynomials)
  - ``T`` -- leading monomials of ``G`` (to discard redundant polynomials)
  - ``strategy`` -- how to sort the critical pairs
  """

  # counters
  cdef int i, j, k
  # signals
  cdef int fails, passes
  # polynomials, clothed and otherwise
  cdef clothed_polynomial cf, cg, ch, cp, cq
  cdef MPolynomial_libsingular f, g, h, p, q, tf, tg, th, tp, tq, tfg, tfh, tpq
  # current critical pair
  cdef tuple pair

  # setup
  cf = G.pop(-1)
  f = cf.value()
  tf = f.lm()
  cdef MPolynomial_libsingular zero = R(0)
  cdef MPolynomial_libsingular one = R(1)

  # some polynomials were removed, so their LTs could be different. fix this
  for pair in P:
    if pair[0].value().parent() != R: pair[0].set_value(R(pair[0].value()))
    if pair[1].value().parent() != R: pair[1].set_value(R(pair[1].value()))

  # STEP 1: eliminate old pairs by Buchberger's lcm criterion
  cdef list Pnew = list()

  for pair in P:

    cp = pair[0]; cq = pair[1]
    p = cp.value(); q = cq.value()

    if q != 0:
      tp = p.lm(); tq = q.lm(); tpq = tp.lcm(tq)

      if (not monomial_divides(tf, tpq)) or tp.lcm(tf) == tpq or tq.lcm(tf) == tpq:
        #print "adding", (tp, tq)
        if strategy=='normal': Pnew.append((cp, cq, lcm_of_critical_pair((cp, cq))))
        else: Pnew.append(pair)

      #else: print (tp,tq,tpq,pair[-1]), "pruned because of", tf

    else:
      if strategy == 'normal': Pnew.append((cp, cq, lcm_of_critical_pair((cp, cq))))
      else: Pnew.append(pair)

  # STEP 2: create list of new pairs
  cdef list C = list()
  cdef list D = list()

  for cg in G:

    if strategy=='sugar': C.append((cf,cg,sug_of_critical_pair((cf,cg))))
    elif strategy=='normal': C.append((cf,cg,lcm_of_critical_pair((cf,cg))))
    elif strategy=='mindeg': C.append((cf,cg,deg_of_critical_pair((cf,cg))))

  # STEP 3: eliminate new pairs by Buchberger's lcm criterion
  i = 0
  while i < len(C):

    pair = C.pop(i)
    #print "considering", pair,
    tfg = lcm_of_critical_pair(pair)

    if tfg == pair[0].lm() * pair[1].lm(): D.append(pair) # relatively prime; postpone removal

    else:

      # first check if unconsidered new critical pairs will prune it
      j = 0
      fails = False
      while j < len(C) and not fails:

        tpq = lcm_of_critical_pair(C[j])

        if monomial_divides(tpq,tfg):
          #print (pair[0].value().lm(), pair[1].value().lm(), lcm_of_critical_pair(pair), pair[-1]), "pruned because of", (C[j][0].value().lm(), C[j][1].value().lm(), lcm_of_critical_pair(C[j]), C[j][-1])
          fails = True

        j += 1

      # now check if considered new critical pairs will prune it
      j = 0
      while j < len(D) and not fails:

        tpq = lcm_of_critical_pair(D[j])

        if monomial_divides(tpq,tfg):
          #print (pair[0].value().lm(), pair[1].value().lm(), lcm_of_critical_pair(pair), pair[-1]), "pruned because of", (D[j][0].value().lm(), D[j][1].value().lm(), lcm_of_critical_pair(D[j]), D[j][-1])
          fails = True

        j += 1

      if not fails:
        D.append(pair)#; print "adding", pair
      #else: print "not added"

  # STEP 4: eliminate new pairs by Buchberger's gcd criterion
  for cg in G:

    g = cg.value(); tg = g.lm()

    if tf.gcd(tg) == one:
      tfg = tf.lcm(tg)
      i = 0

      while i < len(D):

        pair = D[i]; tp = pair[0].value().lm(); tq = pair[1].value().lm()
        tpq = tp.lcm(tq)

        if tpq == tfg:
          D.pop(i)
          #print (tf, tg), "pruned: coprime"

        else: i += 1

  # add new polynomials to basis
  Pnew.extend(D)
  # sort according to strategy
  if strategy == 'sugar': Pnew.sort(key=last_element_then_lcm)
  elif strategy == 'normal': Pnew.sort(key=lcm_of_critical_pair)
  elif strategy == 'mindeg': Pnew.sort(key=deg_of_critical_pair)
  #print [(pair[0].value().lm(), pair[1].value().lm(), lcm_of_critical_pair(pair), pair[-1]) for pair in Pnew]

  # DO NOT REMOVE REDUNDANT ELEMENTS FROM BASIS -- performance suffers

  G.append(cf)
  return Pnew

@cython.profile(True)
def create_order(list w):
  r"""
    Create term ordering acceptable to Singular, using integer weight vector ``w``.
  """
  # first we convert floats to integer
  # this is fine since we are using integer programming to find the ordering
  cdef list wZZ = [ZZ(each) for each in w]
  cdef list M = list()
  M.append(wZZ)

  # now fill in the lower rows of the matrix, to avoid potential ambiguity
  for i in xrange(len(w)-1):
    M.append([1 for k in xrange(i+1,len(w))] + [0 for k in xrange(i+1)])

  return TermOrder(matrix(M))

#TODO fix bug that happens when the number of dynamic iterations is smaller than number of polys
@cython.profile(True)
cpdef tuple dynamic_gb(F, dmax=Infinity, strategy='normal', static=False, minimize_homogeneous=False, \
                       weighted_sugar = 0, use_boundary_vectors=True, use_disjoint_cones=True, \
                       unrestricted=False, random=False, perturbation=False, simplex=False, max_calls=Infinity):
  r"""
    Computes a dynamic Groebner basis of the polynomial ideal generated by ``F``,
    up to degree ``dmax``.

    INPUTS:

      - `F` -- generators of ideal
      - `dmax` -- maximum degree of Groebner basis (for computing `d`-Groebner basis)
      - `strategy` -- one of
        - `mindeg` -- choose critical pairs of minimal degree
        - `normal` -- choose critical pairs of minimal lcm
        - `sugar` -- choose critical pairs of minimal sugar
      - `static` -- whether to use the dynamic algorithm or the static one
      - `minimize_homogeneous` -- whether to keep the weight
        of the homogeneous variable below that of the other variables
      - `weighted_sugar` -- whether the sugar is computed by exponents (0) or by degree (1)
      - `use_boundary_vectors` -- whether to use boundary vectors; setting this and
        `use_disjoint_cones` to `False` is similar to Caboara's old method
      - `use_disjoint_cones` -- whether to use disjoint cones; setting this and
        `use_boundary_vectors` to `False` is similar to Caboara's old method
      - `unrestricted` -- uses Gritzmann-Sturmfels' dynamic algorithm instead of Caboara's
      - `max_calls` -- the maximum number of calls to the dynamic engine
  """
  global sugar_type
  # counters
  cdef int i, j, k
  cdef int calls = 0
  # measures
  cdef int m, n, d, maximum_size_of_intermediate_basis, zero_reductions, number_of_spolynomials
  global rejections, monomials_eliminated, number_of_programs_created, failed_systems
  # signals
  cdef int LTs_changed, sugar_strategy

  # variables related to the polynomials and the basis
  cdef list G = list()
  cdef list generators = list()
  cdef list reducers
  cdef MPolynomial_libsingular p, t
  cdef clothed_polynomial f, g, s, r
  # lists of leading terms
  cdef list LTs = list()
  cdef list newLTs
  # base ring
  cdef MPolynomialRing_libsingular PR

  # variables related to the term ordering
  cdef list current_ordering
  cdef set boundary_vectors
  cdef MixedIntegerLinearProgram lp

  # variables related to the critical pairs
  cdef list P
  cdef tuple Pd

  # variables related to unrestricted dynamic algorithm
  cdef list old_vertices = list()

  # check the strategy first
  if strategy == 'sugar':
    sugar_strategy = True
    sugar_type = weighted_sugar
  else: sugar_strategy = False

  # initialize measures
  rejections = 0; monomials_eliminated = 0; number_of_programs_created = 0; failed_systems = 0
  maximum_size_of_intermediate_basis = 0; zero_reductions = 0; number_of_spolynomials = 0

  # initialize polynomial ring
  PR = F[0].parent()
  n = len(PR.gens())
  cdef clothed_zero = clothed_polynomial(PR(0))
  current_ordering = [1 for k in xrange(n)]

  # set up the linear program and associated variables
  cdef set rejects = set()
  import sage.numerical.backends.glpk_backend as glpk_backend
  lp = new_linear_program()
  # when solving integer programs, perform simplex first, then integer optimization
  # (avoids a GLPK bug IIRC)
  lp.solver_parameter(glpk_backend.glp_simplex_or_intopt, glpk_backend.glp_simplex_then_intopt)

  slp = make_solver()

  # need positive weights
  for k in xrange(n):
    lp.add_constraint(lp[k],min=tolerance_cone)
    #lp.set_min(lp[k],tolerance_cone)
    lp.set_integer(lp[k])
    lp.set_max(lp[k],upper_bound)

  # do we hate the homogenizing variable?
  if minimize_homogeneous:
    for k in xrange(n-1):
      lp.add_constraint(lp[k] - lp[n-1],min=tolerance_cone)

  lp.set_objective(lp.sum([lp[k] for k in xrange(n)]))

  # set up the basis
  m = 0; P = list(); Done = set()
  LTs = list()
  boundary_vectors = None
  reducers = []
  LTs = []

  # clothe the generators
  for p in F:
    f = clothed_polynomial(PR(p))
    generators.append(f)
    if strategy == 'sugar':
      P.append((f,clothed_zero,sug_of_critical_pair((f,clothed_zero))))
    elif strategy == 'normal':
      P.append((f,clothed_zero,lcm_of_critical_pair((f,clothed_zero))))
    elif strategy == 'mindeg':
      P.append((f,clothed_zero,deg_of_critical_pair((f,clothed_zero))))
  m = len(G)

  # initial sort
  if strategy == 'sugar': P.sort(key=last_element_then_lcm)
  elif strategy == 'normal': P.sort(key=lcm_of_critical_pair)
  elif strategy == 'mindeg': P.sort(key=deg_of_critical_pair)

  # main loop
  while len(P) != 0:

    # some diagnostic information
    print "-----------------------"
    #print "current ordering", current_ordering
    print len(P), "critical pairs remaining"
    print len(G), "polynomials in basis"
    maximum_size_of_intermediate_basis = max(maximum_size_of_intermediate_basis, len(G))
    #hp = PR.ideal(LTs).hilbert_polynomial()
    #hs = PR.ideal(LTs).hilbert_series()
    #print "predicted hilbert polynomial", hp
    #print "predicted hilbert series", hs
    #print "predicted hilbert dimension", hp.degree()

    # select critical pairs of minimal degree
    Pd = P.pop(0)

    if d < dmax: # don't go further than requested

      # compute s-polynomials
      #if sugar_strategy: print "current sugar", Pd[len(Pd)-1]
      s = spoly(Pd, generators)
      number_of_spolynomials += 1

      # reduce s-polynomials modulo current basis wrt current order
      if sugar_strategy: r = reduce_polynomial_with_sugar(s, G)
      else: r = reduce_poly(s, reducers)

      # diagnostic
      #print "new polynomial generated",
      #print "leading monomial with current ordering would be", r.value().lm()
      if r.value()==0: zero_reductions += 1

      if r.value() != 0: # add to basis, choose new ordering, update pairs

        #print r
        G.append(r)

        #Start using Perry's restricted algorithm
        if calls == max_calls:
          simplex = False
          unrestricted = False
          perturbation = False
          random = False

        if not static:# and calls < max_calls: # choose a new ordering, coerce to new

          calls += 1

          if unrestricted:
            current_ordering, old_vertices = choose_ordering_unrestricted(G, old_vertices)
          elif random:
            current_ordering = choose_random_ordering(G, current_ordering)
          elif perturbation:
            current_ordering = choose_local_ordering(G, current_ordering)
          elif simplex:
            current_ordering = choose_simplex_ordering(G, current_ordering, slp)
          else:
            current_ordering, lp, boundary_vectors = choose_ordering_restricted(G, LTs[:m], m, current_ordering, \
              lp, rejects, boundary_vectors, use_boundary_vectors, use_disjoint_cones)

          # set up a ring with the current ordering
          #print "current ordering", current_ordering
          PR = PolynomialRing(PR.base_ring(), PR.gens(), order=create_order(current_ordering))
          #print "have ring"
          oldLTs = copy(LTs)
          LTs = list()

          for k in xrange(len(G)): # coerce to new ring
            G[k].set_value(PR(G[k].value()))
            LTs.append(G[k].value().lm())

          if len(oldLTs) > 2 and oldLTs != LTs[:len(LTs)-1]:
            if unrestricted or random or perturbation or simplex:
              changes = [oldLTs[i] == LTs[i] for i in xrange(len(LTs)-1)].count(False)
              print "Changed order:", changes, "changes out of", len(oldLTs), "possible"
              #TODO can I only gm_update stuff that had its LT changed?
              # rebuild P - this is necessary because we changed leading terms
              P = [ Pd for Pd in P if Pd[1].value() == 0 ] #keep unprocessed input polys in queue
              for i in xrange(1, len(G)-1):
                P = gm_update(PR, P, G[:i], LTs[:i], strategy)
            else:
              raise ValueError, "leading terms changed" # this should never happen

          for k in xrange(len(generators)):
            generators[k].set_value(PR(generators[k].value()))

        # setup reduction for next round
        reducers = [G[k].value() for k in xrange(len(G))]
        #print "have new polys and new lts"
        #print PR.term_order()
        P = gm_update(PR, P, G, LTs, strategy)
        #print "updated P"
        m = len(G)


  # done; prune redundant elements
  i = m - 1
  while i > 0:
    p = reducers[i]; t = p.lm()
    m = 0
    while m < i:
      if monomial_divides(t, reducers[m].lm()):
        reducers.pop(m)
        i -= 1
      else: m += 1
    i -= 1

  #The procedure above is not enough to obtain a reduced GB...
  reducers = list(PR.ideal(reducers).interreduced_basis())
  monomials = sum([ len(p.monomials()) for p in reducers ])

  # done; diagnostic and return
  print "final order", current_ordering
  print number_of_spolynomials, "s-polynomials considered (direct count)", len(reducers), "polynomials in basis;", zero_reductions, "zero reductions"
  print "there were at most", maximum_size_of_intermediate_basis, "polynomials in the basis at any one time"
  print "there are", monomials, "monomials in the output basis"
  print rejections, "comparisons were eliminated by rejection"
  print monomials_eliminated, "monomials were eliminated by hypercube vectors"
  print number_of_programs_created, "linear programs were created (including copies)"
  print failed_systems, "failed systems and ", len(rejects), "rejections stored"
  print lp.number_of_constraints(), "constraints are in the linear program"

  #Check that the results are correct
  assert PR.ideal(reducers) == PR.ideal(F), "Output basis generates wrong ideal"
  assert PR.ideal(reducers).gens().is_groebner(), "Output basis is not a GB"
  return reducers, LTs, rejects, G

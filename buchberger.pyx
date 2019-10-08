'''
Base functionality of Buchberger's algorithm for dynamic Buchberger algorithms.

To run, call dynamic_gb.
'''

# cython: profile = False
# cython: boundscheck = False
# cython: wraparound = False
# distutils: language=c++
# distutils: include_dirs=$SAGE_ROOT/local/include/singular
# distutils: libraries=m readline Singular givaro gmpxx gmp

include "types.pxi"
include "polynomials.pxi"
include "stats.pxi"
include "heuristics.pxi"
include "caboara_perry.pxi"
include "unrestricted.pxi"

#Python-level imports

import cython
import time

from copy import copy

from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

@cython.profile(True)
cpdef clothed_polynomial spoly(tuple Pd, int sugar_type):
  """
    Computes an s-polynomial indexed by ``Pd``.

    INPUT:

    - ``Pd`` -- tuple (f,g); we wish to compute spoly(`f,g`) where f and g are
      clothed
  """
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
    #print "building", f.lm()

  #print s
  cs = clothed_polynomial(s, sugar_type); cs.set_sugar(new_sugar)
  return cs

@cython.profile(True)
cpdef clothed_polynomial reduce_polynomial_with_sugar \
  (clothed_polynomial s, list G, int sugar_type):
  """
    Value of ``cf`` changes to reduced polynomial, and new sugar is computed.
    based on algorithm in Cox, Little, and O'Shea, and paper on Sugar by Giovini et al.

    INPUTS:
    - `s` -- s-polynomial to reduce
    - `G` -- list of reducers
  """

  init_time = time.time()
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

  statistics.inc_reduction_time(time.time() - init_time)

  return s

@cython.profile(True)
cpdef clothed_polynomial reduce_poly(clothed_polynomial s, list G):
  """
    Reduces the s-polynomial `s`` modulo the polynomials ``G``
  """
  # polynomials
  cdef MPolynomial_libsingular r

  if s.value() != 0: # no point wasting time
    r = s.value().reduce(G)
    s.set_value(r)

  return s

@cython.profile(True)
cpdef int sug_of_critical_pair(tuple pair, int sugar_type):
  """
    Compute the sugar of a critical pair.
  """
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
      su = sum(u.exponents(as_ETuples=False)[0])
      sv = sum(v.exponents(as_ETuples=False)[0])

    elif sugar_type == 1: # weighted sugar: use degree, determined by ordering
      su = u.degree(); sv = v.degree()

    sfg = max(sf + su, sg + sv)

  return sfg

@cython.profile(True)
cpdef lcm_of_critical_pair(tuple pair):
  """
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
  """
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
cpdef last_element(tuple p):

  return p[len(p)-1] # b/c cython doesn't allow lambda experessions, I think

@cython.profile(True)
cpdef lcm_then_last_element(tuple p):

  return (lcm_of_critical_pair(p), p[len(p)-1])

@cython.profile(True)
cpdef last_element_then_lcm(tuple p):

  return (p[len(p)-1], lcm_of_critical_pair(p))

@cython.profile(True)
cpdef gm_update(MPolynomialRing_libsingular R, list P, list G, list T, \
                strategy, int sugar_type):
  """
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

      if (not monomial_divides(tf, tpq)) or tp.lcm(tf) == tpq \
         or tq.lcm(tf) == tpq:
        #print "adding", (tp, tq)
        if strategy=='normal':

          Pnew.append((cp, cq, lcm_of_critical_pair((cp, cq))))

        else: Pnew.append(pair)

      #else: print (tp,tq,tpq,pair[-1]), "pruned because of", tf

    else:
      if strategy == 'normal':

        Pnew.append((cp, cq, lcm_of_critical_pair((cp, cq))))

      else: Pnew.append(pair)

  # STEP 2: create list of new pairs
  cdef list C = list()
  cdef list D = list()

  for cg in G:

    if strategy=='sugar':

      C.append((cf,cg,sug_of_critical_pair((cf,cg), sugar_type)))

    elif strategy=='normal': C.append((cf,cg,lcm_of_critical_pair((cf,cg))))
    elif strategy=='mindeg': C.append((cf,cg,deg_of_critical_pair((cf,cg))))

  # STEP 3: eliminate new pairs by Buchberger's lcm criterion
  i = 0
  while i < len(C):

    pair = C.pop(i)
    #print "considering", pair,
    tfg = lcm_of_critical_pair(pair)

    if tfg == pair[0].lm() * pair[1].lm():

      D.append(pair) # relatively prime; postpone removal

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

cpdef list rebuild_queue(list G, list LMs, list P, str strategy, int sugar_type):

  cdef list Pnew = [ Pd for Pd in P if Pd[1].value() == 0 ]
  cdef int i
  cdef MPolynomialRing_libsingular R = G[0].value().parent()
  for i in xrange(1, 1+len(G)):
    Pnew = gm_update(R, Pnew, G[:i], LMs[:i], strategy, sugar_type)

  return Pnew

#List of algorithms:
# static, caboara, caboara-perry, gritzmann-sturmfels, random, perturbation,
# simplex, regrets.
@cython.profile(True)
cpdef tuple dynamic_gb \
  (F, dmax=Infinity, strategy='sugar', \
   minimize_homogeneous=False, weighted_sugar = 0, algorithm='static', \
   max_calls=Infinity, itmax=Infinity, print_results=False, \
   print_candidates = False, heuristic='hilbert', initial_ordering='grevlex',
   dynamic_period = 10, reducer='classical'):
  """
    Computes a dynamic Groebner basis of the polynomial ideal generated by
      ``F``, up to degree ``dmax``.

    INPUTS:

      - `F` -- generators of ideal
      - `algorithm` -- one of [ static, cabaora, caboara-perry, gritzmann-sturmfels, random, perturbation, simplex, regrets ]
      - `dmax` -- maximum degree of Groebner basis (for computing `d`-Groebner
        basis)
      - `strategy` -- one of
        - `mindeg` -- choose critical pairs of minimal degree
        - `normal` -- choose critical pairs of minimal lcm
        - `sugar` -- choose critical pairs of minimal sugar
      - `minimize_homogeneous` -- whether to keep the weight
        of the homogeneous variable below that of the other variables
      - `weighted_sugar` -- whether the sugar is computed by exponents (0) or by
        degree (1)
      - `max_calls` -- the maximum number of calls to the dynamic engine
      - `itmax` -- run for `itmax` iterations only and return ordering
      - `print_results` -- print GB data (size, number of monomials, degree, etc)
      - `print_candidates` -- print candidate LTs (for debugging)
      - `heuristic` -- heuristic function to use. hilbert, betti, or mixed
      - `reducer` -- either 'classical' (Buchberger algorithm) or 'F4'
  """
  statistics.set_print_results(print_results)
  statistics.reset_all_stats()
  cdef int sugar_type
  # counters
  cdef int i, j, k
  cdef int calls = 0
  cdef int iteration_count = 0
  # measures
  cdef int m, n, d = 0
  # signals
  cdef int LTs_changed, sugar_strategy

  # variables related to the polynomials and the basis
  cdef list G = list()
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

  # check the strategy first
  if strategy == 'sugar':
    sugar_strategy = True
    sugar_type = weighted_sugar
  else: sugar_strategy = False

  # initialize polynomial ring
  PR = F[0].parent()
  n = len(PR.gens())
  if initial_ordering == 'random':
      current_ordering = [ randint(1, 10) for k in xrange(n) ]
  elif initial_ordering == 'grevlex':
      current_ordering = [1 for k in xrange(n)]
  else:
      current_ordering = initial_ordering
  PR = PolynomialRing(PR.base_ring(), PR.gens(), order=create_order(current_ordering))
  cdef clothed_zero = clothed_polynomial(PR(0), sugar_type)


  # set up the linear program and associated variables
  cdef set rejects = set()
  lp = new_linear_program(n = n)

  #State for additional algorithms
  if algorithm == 'caboara':
    use_disjoint_cones = False
    use_boundary_vectors = False
  if algorithm == 'caboara-perry':
    use_disjoint_cones = True
    use_boundary_vectors = True
  if algorithm == 'gritzmann-sturmfels':
    old_polyhedron = Polyhedron(rays=(-identity_matrix(n)).rows())
  if algorithm == 'simplex':
    slp = make_solver(n)
    init_linear_program(slp, n)
  cdef list constraints = []
  cdef list vertices = []

  cdef int prev_hilbert_degree = -1

  # set up the basis
  m = 0; P = list(); Done = set()
  LTs = list()
  boundary_vectors = None
  reducers = []
  LTs = []

  # clothe the generators
  for p in F:
    f = clothed_polynomial(PR(p), sugar_type)
    if strategy == 'sugar':
      P.append((f,clothed_zero,\
                sug_of_critical_pair((f,clothed_zero), sugar_type)))
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
    #print "-----------------------"
    #print "current ordering", current_ordering
    #print len(P), "critical pairs remaining"
    #print len(G), "polynomials in basis"
    #maximum_size_of_intermediate_basis = max(maximum_size_of_intermediate_basis, len(G))
    statistics.update_maximum_intermediate_basis(len(G))
    #hp = PR.ideal(LTs).hilbert_polynomial()
    #hs = PR.ideal(LTs).hilbert_series()
    #print "predicted hilbert polynomial", hp
    #print "predicted hilbert series", hs
    #print "predicted hilbert dimension", hp.degree()

    # select critical pairs of minimal sugar / degree / ...
    Pd = P.pop(0)

    #Stop here and return ordering if asked
    if iteration_count >= itmax:
      return (current_ordering, )
    else:
      iteration_count += 1

    if d < dmax: # don't go further than requested

      # compute s-polynomials
      #if sugar_strategy: print "current suga", Pd[len(Pd)-1]
      s = spoly(Pd, sugar_type)
      statistics.inc_spolynomials()

      # reduce s-polynomials modulo current basis wrt current order
      if sugar_strategy: r = reduce_polynomial_with_sugar(s, G, sugar_type)
      else: r = reduce_poly(s, reducers)

      # diagnostic
      #print "new polynomial generated",
      #print "leading monomial with current ordering would be", r.value().lm()
      if r.value()==0: statistics.inc_zero_reductions()

      if r.value() != 0: # add to basis, choose new ordering, update pairs

        G.append(r)

        #Start using Perry's restricted algorithm when limit to calls of
        #unrestricted algorithm has been reached
        if calls == max_calls:
          algorithm = 'caboara-perry'

        #TODO decide whether to take max_calls or not!
        #TODO maybe use some other criterion to activate dynamic algorithms
        if algorithm != 'static':# and calls < max_calls: # choose a new ordering, coerce to new

          calls += 1
          dynamic_time = time.time()

          if algorithm == 'gritzmann-sturmfels':
            current_ordering, old_polyhedron, prev_hilbert_degree = \
              choose_ordering_unrestricted(G, old_polyhedron, heuristic, \
                                           len(P), prev_hilbert_degree)
          elif algorithm == 'random':
            current_ordering, prev_hilbert_degree = choose_random_ordering \
                (G, current_ordering, heuristic, len(P), prev_hilbert_degree)
          elif algorithm == 'perturbation':
            if iteration_count % dynamic_period == 1:
              current_ordering, prev_hilbert_degree = choose_perturbation_ordering \
                (G, current_ordering, heuristic, len(P), prev_hilbert_degree)
          elif algorithm == 'simplex':
            current_ordering, vertices, prev_hilbert_degree = \
              choose_simplex_ordering(G, current_ordering, slp, vertices, \
                                      heuristic, len(P), prev_hilbert_degree)
          elif algorithm == 'regrets':
            current_ordering, lp, boundary_vectors, constraints, changed = \
              choose_regrets_ordering(G, current_ordering, constraints, lp, \
                                   heuristic)
          else:
            current_ordering, lp, boundary_vectors = \
              choose_ordering_restricted(G, LTs[:m], m, current_ordering, \
                                         lp, rejects, boundary_vectors,
                                         use_boundary_vectors,
                                         use_disjoint_cones, print_candidates, \
                                         heuristic)

          statistics.inc_dynamic_overhead(time.time() - dynamic_time)

          # set up a ring with the current ordering
          #print "current ordering", current_ordering
          PR = PolynomialRing(PR.base_ring(), PR.gens(), order=create_order(current_ordering))
          #print "have ring"
          oldLTs = copy(LTs)
          LTs = list()

          for k in xrange(len(G)): # coerce to new ring
            G[k].set_value(PR(G[k].value()))
            LTs.append(G[k].value().lm())

          if len(oldLTs) > 0 and oldLTs != LTs[:len(LTs)-1]:
            if algorithm in [ 'gritzmann-sturmfels', 'random', 'perturbation', 'simplex' ]:
              queue_time = time.time()
              if iteration_count % dynamic_period == 1:
                P = rebuild_queue(G[:len(G)-1], LTs[:len(LTs)-1], P, strategy, sugar_type)
              statistics.inc_queue_time(time.time() - queue_time)
            elif algorithm == 'regrets':
              queue_time = time.time()
              P = gm_update(PR, P, G[:len(G)-1], LTs[:len(LTs)-1], strategy, \
                            sugar_type) #do update w.r.t new poly
              statistics.inc_queue_time(time.time() - queue_time)
            else:
              raise ValueError, "leading terms changed" # this should never happen

        # setup reduction for next round
        reducers = [G[k].value() for k in xrange(len(G))]
        #print "have new polys and new lts"
        #print PR.term_order()
        queue_time = time.time()
        P = gm_update(PR, P, G, LTs, strategy, sugar_type)
        statistics.inc_queue_time(time.time() - queue_time)
        #print "queue", len(P)
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

  #Update and print statistics on the execution
  statistics.update_running_time()
  statistics.set_number_of_rejects(len(rejects))
  statistics.set_number_of_constraints(lp.number_of_constraints())
  statistics.update_basis_data(reducers)
  statistics.brief_report()

  #print current_ordering

  #Check that the results are correct
  #assert PR.ideal(reducers).gens().is_groebner(), "Output basis is not a GB"
  #assert PR.ideal(reducers) == PR.ideal(F), "Output basis generates wrong ideal"
  return reducers, LTs, rejects, G

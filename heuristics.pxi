cpdef int hs_heuristic(tuple f, tuple g):
  r"""
    Implements the Hilbert heuristic recommended by Caboara in his 1993 paper.
    It first compares the degrees of the Hilbert polynomials;
    if one polynomial's degree is smaller than the other,
    then it returns the difference b/w the degree for f's hp and g's.
    Otherwise, it compares the Hilbert series of the polynomials;
    if the two are equal, it returns 0; otherwise,
    it returns the trailing coefficient of the numerator of the difference b/w
    f's hs and g's. See Caboara's paper for details
  """

  if f[0].degree() == g[0].degree():

    if f[1] == g[1]: return 0
    else:
      C = (f[1]-g[1]).numerator().coefficients()
      return C[len(C)-1]

  else: return f[0].degree() - g[0].degree()

@cython.profile(True)
cpdef list sort_CLTs_by_Hilbert_heuristic(MPolynomialRing_libsingular R, \
                                          list current_Ts, list CLTs):
  r"""
    Sorts the compatible leading monomials using the Hilbert heuristic
    recommended by Caboara in his 1993 paper.
    Preferable leading monomials appear earlier in the list.

    INPUTS:
    - `R` -- current ring (updated according to latest ordering!)
    - `current_Ts` -- current list of leading terms
    - `CLTs` -- compatible leading terms from which we want to select leading
      term of new polynomial
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
           tup) for tup in CLTs]
  CLTs.sort(cmp=hs_heuristic) # sort according to hilbert heuristic
  #print CLTs

  return CLTs

cpdef list min_CLT_by_Hilbert_heuristic(MPolynomialRing_libsingular R, \
                                        list CLTs):
  r"""
    Sorts the compatible leading monomials using the Hilbert heuristic
    recommended by Caboara in his 1993 paper.
    Preferable leading monomials appear earlier in the list.

    INPUTS:
    - `R` -- current ring (updated according to latest ordering!)
    - `CLTs` -- compatible leading terms from which we want to select
      leading term of new polynomial

    OUTPUTS:
    - the preferred list of leading monomials from CLTs
  """

  cdef list leads #lists of leading monomials in CLTs

  CLTs = [(R.ideal(leads).hilbert_polynomial(),
           R.ideal(leads).hilbert_series(),
           leads) for leads in CLTs]
  CLTs.sort(cmp=hs_heuristic)
  #TODO this could be optimized, we are sorting just to find the minimum...
  return CLTs[0][2]

cpdef list min_weights_by_Hilbert_heuristic(MPolynomialRing_libsingular R, \
                                            list CLTs):
  r"""
    Sorts the compatible leading monomials using the Hilbert heuristic
    recommended by Caboara in his 1993 paper.
    Preferable leading monomials appear earlier in the list.

    INPUTS:
    - `R` -- current ring (updated according to latest ordering!)
    - `current_Ts` -- current list of leading terms
    - `CLTs` -- compatible leading terms from which we want to select leading
      term of new polynomial

    OUTPUTS:
    - the preferred weight vector from CLTs
  """

  cdef tuple leads #lists of leading monomials in CLTs

  CLTs = [(R.ideal(leads[1]).hilbert_polynomial(),
           R.ideal(leads[1]).hilbert_series(),
           leads[0]) for leads in CLTs]
  CLTs.sort(cmp=hs_heuristic)

  #TODO this could be optimized, we are sorting just to find the minimum...
  return CLTs[0][2]

##Implementation of Betti-based heuristics

cpdef bool is_edge(int i, int j, list LMs):
  r"""
  Returns true iff (i, j) is an edge in the Buchberger graph of LMs
  """
  cdef MPolynomialRing_libsingular R = LMs[0].parent()
  cdef int k, m
  for k in xrange(len(LMs)):
    if k == i or k == j:
      continue
    for v in R.gens():
      #TODO careful here, this is not quite it!
      m = max([LMs[i].degree(v), LMs[j].degree(v)])
      if LMs[k].degree(v) >= m:
        if m > 0 or LMs[k].degree(v) > m:
          break #LMs[k] does not strictly divide in every variable lcm(LMs[i], LMs[j])
    else: #LMs[k] strictly divides lcm(LMs[i], LMs[j]) in every variable
      return False
  return True

cpdef int graph_edges(list LMs):
  cdef int num_edges = 0
  cdef int j, i
  for j in xrange(len(LMs)):
    for i in xrange(j):
      if is_edge(i, j, LMs):
        num_edges += 1
  return num_edges

cpdef int betti_heuristic(tuple f, tuple g):
  return f[0] - g[0]

cpdef int hilbert_betti_heuristic(tuple f, tuple g):

  if f[0].degree() == g[0].degree():
    # Break Hilbert degree ties by Betti number approximation
    return f[1] - g[1]

  return f[0].degree() - g[0].degree()

cpdef list sort_CLTs_by_heuristic(list CLTs, str heuristic, bool use_weights, \
                                  int prev_betti=-1, int prev_hilb=-1):

  # if use_weights, CLTs is a list of tuples (candidate lts, weight vector)
  # otherwise, it is a list of candidate lts
  cdef list L, old_order
  cdef MPolynomialRing_libsingular R = CLTs[0][0][0].parent()
  if heuristic == 'hilbert':

    if use_weights:
      L = [ (R.ideal(LTs[0]).hilbert_polynomial(),
             R.ideal(LTs[0]).hilbert_series(),
             LTs) for LTs in CLTs ]
    else:
      L = [ (R.ideal(LTs).hilbert_polynomial(),
             R.ideal(LTs).hilbert_series(),
             LTs) for LTs in CLTs ]
    L.sort(cmp=hs_heuristic)
    return L

  elif heuristic == 'betti':

    old_order = [ ((), (), CLTs[0]) ]
    if use_weights:
      L = [ (graph_edges(LTs[0]), (), LTs) for LTs in CLTs ]
    else:
      L = [ (graph_edges(LTs), (), LTs) for LTs in CLTs ]
    L.sort(cmp=betti_heuristic)
    if prev_betti >= 0 and prev_betti < L[0][0]:
      return old_order
    return L

  elif heuristic == 'mixed':

    old_order = [ ((), (), CLTs[0]) ]
    if use_weights:
      L = [ (R.ideal(LTs[0]).hilbert_polynomial(),
             graph_edges(LTs[0]),
             LTs) for LTs in CLTs ]
    else:
      L = [ (R.ideal(LTs).hilbert_polynomial(),
             graph_edges(LTs),
             LTs) for LTs in CLTs ]
    L.sort(cmp=hilbert_betti_heuristic)
    if prev_hilb >= L[0][0] and prev_betti >= 0 and prev_betti < L[0][1]:
      return old_order
    return L

  raise ValueError("Invalid heuristic function: " + heuristic)

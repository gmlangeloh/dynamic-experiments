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
  cdef list CLTs = [(R.ideal(current_Ts + [tup[1]]).hilbert_polynomial(),
                     R.ideal(current_Ts + [tup[1]]).hilbert_series(),
                     tup)
                    for tup in CLTs]
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

  cdef list CLTs = [(R.ideal(leads).hilbert_polynomial(),
                     R.ideal(leads).hilbert_series(),
                     leads)
                    for leads in CLTs]
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

  cdef list CLTs = [(R.ideal(leads[1]).hilbert_polynomial(),
                     R.ideal(leads[1]).hilbert_series(),
                     leads[0])
                    for leads in CLTs]
  CLTs.sort(cmp=hs_heuristic)

  #TODO this could be optimized, we are sorting just to find the minimum...
  return CLTs[0][2]

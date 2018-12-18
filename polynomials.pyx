from sage.libs.singular.decl cimport p_DivisibleBy

cdef class clothed_polynomial:

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

  cpdef MPolynomial_libsingular value(self): return self.f

  cpdef set_value(self, MPolynomial_libsingular f): self.f = f

  cpdef lm(self): return self.f.lm()

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

from libcpp cimport bool
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular


#Defining fast monomial operations with fewer checks
cpdef int monomial_divides(MPolynomial_libsingular t, MPolynomial_libsingular u)

cpdef int indivisible(tuple t, tuple u)

#Polynomial class used throughout the dynamic algorithms
cdef class clothed_polynomial:

  r"""
  We use clothed polynomials to store information about polynomials,
  and to make it easier to change the ordering.
  """

  # internal data

  cdef MPolynomial_libsingular f # the polynomial that we have clothed
  cdef int sugar # the sugar of this polynomial, as computed by us

  cdef is_equal(self, clothed_polynomial other)

  # methods related to the polynomial

  cpdef MPolynomial_libsingular value(self)

  cpdef set_value(self, MPolynomial_libsingular f)

  cpdef lm(self)

  cpdef set_sugar(self, int s)

  cpdef int get_sugar(self)

from types cimport *

cdef class Stats:

  cdef bool print_results

  #Statistics
  cdef int rejections
  cdef int monomials_eliminated
  cdef int number_of_programs_created
  cdef int failed_systems
  cdef int maximum_size_of_intermediate_basis
  cdef int zero_reductions
  cdef int number_of_spolynomials

  #Updating statistics
  cpdef void reset_all_stats(self)
  cpdef void inc_rejections(self)
  cpdef void inc_monomials_eliminated(self)
  cpdef void inc_programs_created(self)
  cpdef void inc_failed_systems(self)
  cpdef void inc_maximum_intermediate_basis(self)
  cpdef void inc_zero_reductions(self)
  cpdef void inc_spolynomials(self)

  #Various methods
  cpdef void set_print_results(self, bool value)
  cpdef void report(self)

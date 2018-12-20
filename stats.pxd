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
  cdef int number_of_rejects
  cdef int number_of_constraints

  #Solution data
  cdef list ordering
  cdef int basis_size
  cdef int basis_monomials

  #Updating statistics
  cpdef void reset_all_stats(self)
  cpdef void inc_rejections(self)
  cpdef void inc_monomials_eliminated(self)
  cpdef void inc_programs_created(self)
  cpdef void inc_failed_systems(self)
  cpdef void update_maximum_intermediate_basis(self, int value)
  cpdef void inc_zero_reductions(self)
  cpdef void inc_spolynomials(self)
  cpdef void set_number_of_rejects(self, int value)
  cpdef void set_number_of_constraints(self, int value)

  cpdef void update_basis_data(self, list basis)

  #Various methods
  cpdef void set_print_results(self, bool value)
  cpdef void report(self)

cdef Stats statistics

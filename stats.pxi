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

  cdef float dynamic_overhead
  cdef float running_time

  #Solution data
  cdef list ordering
  cdef int basis_size
  cdef int basis_monomials
  cdef int basis_max_degree

  def __cinit__(self):

    self.print_results = True

  def __init__(self):

    self.reset_all_stats()

  cpdef void set_print_results(self, bool value):

    self.print_results = value

  cpdef void reset_all_stats(self):

    self.rejections = 0
    self.monomials_eliminated = 0
    self.number_of_programs_created = 0
    self.failed_systems = 0
    self.maximum_size_of_intermediate_basis = 0
    self.zero_reductions = 0
    self.number_of_spolynomials = 0
    self.number_of_rejects = 0
    self.number_of_constraints = 0
    self.dynamic_overhead = 0.0
    self.running_time = 0.0

    self.ordering = []
    self.basis_size = 0
    self.basis_monomials = 0
    self.basis_max_degree = 0

  cpdef void inc_dynamic_overhead(self, float val): self.dynamic_overhead += val

  cpdef void update_running_time(self, float time): self.running_time = time

  cpdef void inc_rejections(self): self.rejections += 1

  cpdef void inc_monomials_eliminated(self): self.monomials_eliminated += 1

  cpdef void inc_programs_created(self): self.number_of_programs_created += 1

  cpdef void inc_failed_systems(self): self.failed_systems += 1

  cpdef void update_maximum_intermediate_basis(self, int value):
    self.maximum_size_of_intermediate_basis = \
        max(self.maximum_size_of_intermediate_basis, value)

  cpdef void inc_zero_reductions(self): self.zero_reductions += 1

  cpdef void inc_spolynomials(self): self.number_of_spolynomials += 1

  cpdef void set_number_of_rejects(self, int value):
    self.number_of_rejects = value

  cpdef void set_number_of_constraints(self, int value):
    self.number_of_constraints = value

  cpdef void update_basis_data(self, list basis):
    self.basis_size = len(basis)
    self.basis_monomials = sum([ len(p.monomials()) for p in basis ])
    self.basis_max_degree = max([ p.total_degree(True) for p in basis ])

    if len(basis) > 0:
        self.ordering = list(basis[0].parent().term_order().weights())

  cpdef void report(self):
    if self.print_results:
      print("final order", self.ordering)
      print(self.number_of_spolynomials, "s-polynomials considered (direct count)", self.basis_size, "polynomials in basis;", self.zero_reductions, "zero reductions")
      print("there were at most", self.maximum_size_of_intermediate_basis, "polynomials in the basis at any one time")
      print("there are", self.basis_monomials, "monomials in the output basis")
      print(self.rejections, "comparisons were eliminated by rejection")
      print(self.monomials_eliminated, "monomials were eliminated by hypercube vectors")
      print(self.number_of_programs_created, "linear programs were created (including copies)")
      print(self.failed_systems, "failed systems and ", self.number_of_rejects, "rejections stored")
      print(self.number_of_constraints, "constraints are in the linear program")

  cpdef void brief_report(self):
    if self.print_results:
        print(self.running_time, self.dynamic_overhead, self.basis_size, \
            self.basis_monomials, self.basis_max_degree, \
            self.number_of_spolynomials)

#I will use this stats instance everywhere
statistics = Stats()

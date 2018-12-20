cdef class Stats:

  def __init__(self):

    self.print_results = True
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

  cpdef void inc_rejections(self): self.rejections += 1

  cpdef void inc_monomials_eliminated(self): self.monomials_eliminated += 1

  cpdef void inc_programs_created(self): self.number_of_programs_created += 1

  cpdef void inc_failed_systems(self): self.failed_systems += 1

  cpdef void update_maximum_intermediate_basis(self, int value):
    self.maximum_size_of_intermediate_basis = \
        max(self.maximum_size_of_intermediate_basis, value)

  cpdef void inc_zero_reductions(self): self.zero_reductions += 1

  cpdef void inc_spolynomials(self): self.number_of_spolynomials += 1

  cpdef void report(self):

    if self.print_results:
      #TODO this is obviously missing stuff

      print "final order", current_ordering
      print self.number_of_spolynomials, "s-polynomials considered (direct count)", len(reducers), "polynomials in basis;", self.zero_reductions, "zero reductions"
      print "there were at most", self.maximum_size_of_intermediate_basis, "polynomials in the basis at any one time"
      print "there are", monomials, "monomials in the output basis"
      print self.rejections, "comparisons were eliminated by rejection"
      print self.monomials_eliminated, "monomials were eliminated by hypercube vectors"
      print self.number_of_programs_created, "linear programs were created (including copies)"
      print self.failed_systems, "failed systems and ", len(rejects), "rejections stored"
      print lp.number_of_constraints(), "constraints are in the linear program"
      print "Reinserted", change_count, "times with changes."
      print "Reinserted", no_change_count, "times with no changes."

#I will use this stats instance everywhere
statistics = Stats()

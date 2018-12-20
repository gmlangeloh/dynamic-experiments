from types cimport *
from stats cimport statistics
from polynomials cimport monomial_divides
from heuristics cimport sort_CLTs_by_Hilbert_heuristic

cpdef tuple choose_ordering_restricted \
    (list G, list current_Ts, int mold, list current_ordering, \
     MixedIntegerLinearProgram lp, set rejects, set bvs, int use_bvs, \
     int use_dcs, bool print_candidates)

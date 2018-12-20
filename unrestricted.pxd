from types cimport *
from stats cimport statistics
from polynomials cimport clothed_polynomial, monomial_divides, create_order
from caboara_perry cimport choose_ordering_restricted, new_linear_program, \
    tolerance_cone
from heuristics cimport hs_heuristic, min_CLT_by_Hilbert_heuristic, \
    min_weights_by_Hilbert_heuristic

cpdef GLPKBackend make_solver(int n)

cpdef tuple choose_ordering_unrestricted(list G, list old_vertices)

cpdef tuple choose_simplex_ordering(list G, list current_ordering, \
                                    GLPKBackend lp, list vertices, \
                                    int iterations = *)

cpdef list choose_random_ordering(list G, list current_ordering, \
                                  int iterations = *)

cpdef list choose_local_ordering(list G, list current_ordering, \
                                 int iterations = *)

cpdef tuple choose_cone_ordering \
    (list G, list current_ordering, list constraints, \
     MixedIntegerLinearProgram lp, set rejects, set bvs, int use_bvs, \
     int use_dcs)

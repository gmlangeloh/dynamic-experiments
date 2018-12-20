from types cimport *

cpdef int hs_heuristic(tuple f, tuple g)

cpdef list sort_CLTs_by_Hilbert_heuristic(MPolynomialRing_libsingular R, \
                                          list current_Ts, list CLTs)

cpdef list min_CLT_by_Hilbert_heuristic(MPolynomialRing_libsingular R, \
                                        list CLTs)

cpdef list min_weights_by_Hilbert_heuristic(MPolynomialRing_libsingular R, \
                                            list CLTs)

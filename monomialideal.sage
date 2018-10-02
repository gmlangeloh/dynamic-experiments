r'''
Defines a comparable MonomialIdeal appearing in dynamic algorithms.
Uses Hilbert Function as comparison heuristic.
'''

class MonomialIdeal:
    r'''
    A MonomialIdeal is used to keep an ideal with its associated Hilbert
    Function. MonomialIdeals can be compared with respect to their Hilbert
    Polynomials.
    '''
    def __init__(self, I, w = None, vertices = None):
        self._ideal = I
        self._ring = I.ring()
        if w is not None:
            self._weights = vector(w)
        else:
            self._weights = None
        self._n = I.ring().ngens()
        self._vertices = vertices
        self._hilbert = None

    def hilbert(self):
        if self._hilbert is None:
            self._hilbert = self._ideal.hilbert_polynomial()
        return self._hilbert

    def vertices(self):
        return self._vertices

    def weights(self):
        return self._weights

    def __cmp__(self, other):
        r'''
        self < other iff the degree of its Hilbert Polynomial is smaller or
        degrees are the same and leading coefficient is smaller.
        '''
        HP1 = self.hilbert()
        HP2 = other.hilbert()
        if HP1.degree() < HP2.degree():
            return -1
        elif HP1.degree() == HP2.degree():
            if HP1.lc() < HP2.lc():
                return -1
            elif HP1.lc() == HP2.lc():
                return 0
        return 1

    def _lm_from_weights(self, g):
        max_val = None
        max_mon = None
        for v in g.graph().vertices():
            val = vector(v) * self.weights()
            if max_val is None or val > max_val:
                max_val = val
                max_mon = v
        return prod([ self._ring.gens()[i]^max_mon[i] \
                      for i in range(self._n) ]), max_mon

    def update(self, g):
        old_gens = self._ideal.gens()
        new_gen, new_v = self._lm_from_weights(g)
        self._ideal = ideal(old_gens + [ new_gen ])
        self._vertices += [ new_v ]
        self._hilbert = None

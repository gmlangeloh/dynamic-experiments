'''
Basic implementation of F4 reduction in Sage.
It is not meant to be very efficient, but hopefully can compete with
Caboara and Perry's (2014) implementation of classical reductions.
'''

cpdef symbolic_preprocessing (list L, list G):
  pass

cpdef polynomials_from_matrix (M):
  pass

cpdef reduce_F4 (list L, list G):
  '''
  Build a matrix from a list of pairs L and reduces w.r.t. G.
  For now, I am assuming the normal selection strategy - this is relevant
  here because it means I don't have to worry about updating sugars.
  Faugère reports that the normal strategy works better in his original F4 paper.
  '''
  M = symbolic_preprocessing(L, G) #Build the matrix
  Mred = M.rref() #Do row reduction
  return polynomials_from_matrix(M)

cpdef select_pairs_normal_F4 (list P):
  '''
  Select pairs to be reduced in F4 using the normal strategy, that is,
  picks the critical pairs with minimal degree.
  '''
  if not P: #P is already empty
    return

  cdef int min_deg = P[0][2] #Degree of the first element of the S-polynomial
  #queue. Note that this works because we assume P is already sorted.

  cdef int i = 0
  cdef list L = []
  cdef MPolynomial_libsingular f, g, tf, tg, tfg, s1, s2
  cdef MPolynomialRing_libsingular R
  while P[i][2] == min_deg:

    #compute both "branches" of the S-polynomial and add to the reducer list
    f = P[i][0]
    g = P[i][1]
    tf = f.lm()
    tg = g.lm()
    tgf = tf.lcm(tg)
    R = tfg.parent()
    s1 = R.monomial_quotient(tfg, tf) * f
    s2 = R.monomial_quotient(tfg, tg) * g
    #TODO for efficiency, compute the list of terms here. These are useful in
    #symbolic preprocessing

    L.append((s1, s2))
    i += 1
  return L

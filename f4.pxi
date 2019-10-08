'''
Basic implementation of F4 reduction in Sage.
It is not meant to be very efficient, but hopefully can compete with
Caboara and Perry's (2014) implementation of classical reductions.
'''

from sage.matrix.matrix_modn_sparse cimport Matrix_modn_sparse

cpdef tuple symbolic_preprocessing (list L, set todo, list G):

  #Step 1: Finish computing L, the list of polynomials that will appear as rows
  #of the matrix
  cdef MPolynomial_libsingular f, g #Polynomials
  cdef MPolynomial_libsingular m, n #Monomials
  cdef MPolynomialRing_libsingular R = G[0].parent()
  while todo:
    m = todo.pop()
    #Check if m is top-reducible by G. If it is, add corresponding reducer
    for g in G:
      tg = g.lm()
      if monomial_divides(tg, m):
        n = R.monomial_quotient(m, tg)
        f = n * g
        L.append(f)
        todo.update(f.monomials())
        break

  #Step 2: Build the F4 matrix from L
  #TODO this is inefficient! Do this better later.

  #find set of all monomials in L, sort them (decreasing order)
  cdef set monomial_set = set()
  for f in L:
    monomial_set.update(f.monomials())
  cdef list monomial_list = list(monomial_set)
  monomial_list.sort(reverse=True) #Sort in descending order

  #then make the relation indices <--> coefs, create dictionary
  cdef dict indices_to_coefs = {}
  cdef int i, j
  for i in range(len(F)):
    g = F[i]
    for m in g.monomials():
      j = monomial_list.index(m)
      indices_to_coefs[(i, j)] = g.monomial_coefficient(m)

  #and finally write this dict as a sparse Sage matrix
  cdef Matrix_modn_sparse M = matrix(R.base_ring(), len(F), len(monomial_list),
                                     indices_to_coefs, sparse=True)

  return M, monomial_list

cpdef list polynomials_from_matrix (Matrix_modn_sparse M,
                                    list monomial_list,
                                    MPolynomialRing_libsingular R):

  cdef list new_polys = [ R(0) ] * M.nrows()
  cdef int i, j
  for i, j in M.dict():
    new_polys[i] += M[i][j] * monomial_list[j]

  return new_polys

cpdef list reduce_F4 (list L, set todo, list G):
  '''
  Build a matrix from a list of pairs L and reduces w.r.t. G.
  For now, I am assuming the normal selection strategy - this is relevant
  here because it means I don't have to worry about updating sugars.
  Faugère reports that the normal strategy works better in his original F4 paper.
  '''
  M, monomials = symbolic_preprocessing(L, todo, G) #Build the matrix
  Mred = M.rref() #Do row reduction
  return polynomials_from_matrix(M, monomials, G[0].parent())

cpdef tuple select_pairs_normal_F4 (list P):
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
  cdef set todo = set()
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

    #taking reducers from these polynomials
    #slicing off the first monomial of their lists (i.e., their lms)
    todo.update(s1.monomials()[1:], s2.monomials()[1:])

    L.append((s1, s2))
    i += 1
  return L, todo

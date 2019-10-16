'''
Basic implementation of F4 reduction in Sage.
It is not meant to be very efficient, but hopefully can compete with
Caboara and Perry's (2014) implementation of classical reductions.
'''

from sage.matrix.matrix_modn_sparse cimport Matrix_modn_sparse

@cython.profile(True)
cpdef tuple symbolic_preprocessing (list L, set todo, list G):

  #Step 1: Finish computing L, the list of polynomials that will appear as rows
  #of the matrix
  cdef clothed_polynomial sg
  cdef MPolynomial_libsingular f, g #Polynomials
  cdef MPolynomial_libsingular tg, m, n #Monomials
  cdef MPolynomialRing_libsingular R = L[0].parent()
  cdef set done = set([ f.lm() for f in L ])
  while todo:
    m = todo.pop()
    if m in done:
      continue
    done.add(m)
    #Check if m is top-reducible by G. If it is, add corresponding reducer
    for sg in G:
      g = sg.value()
      tg = g.lm()
      if monomial_divides(tg, m):
        n = R.monomial_quotient(m, tg)
        f = n * g
        L.append(f)
        todo.update([ n for n in f.monomials() if n not in done ])
        break

  #Step 2: Build the F4 matrix from L

  #find set of all monomials in L, sort them (decreasing order)
  cdef set monomial_set = set()
  for f in L:
    monomial_set.update(f.monomials())
  cdef list monomial_list = list(monomial_set)
  monomial_list.sort(reverse=True) #Sort in descending order

  cdef dict monomial_indices = {}
  cdef int i, j
  for i in range(len(monomial_list)):
    monomial_indices[monomial_list[i]] = i

  #TODO building dictionary of indices is still relatively expensive.
  cdef dict indices_to_coefs = {}
  for i in range(len(L)):
    g = L[i]
    for m in g.monomials():
      j = monomial_indices[m]
      indices_to_coefs[(i, j)] = g.monomial_coefficient(m)

  #and finally write this dict as a sparse Sage matrix
  cdef Matrix_modn_sparse M = matrix(R.base_ring(), len(L), len(monomial_list),
                                     indices_to_coefs, sparse=True)

  return M, L, monomial_list

@cython.profile(True)
cpdef add_polynomials_from_matrix (Matrix_modn_sparse M,
                                   list reducers,
                                   list monomial_list,
                                   MPolynomialRing_libsingular R,
                                   list G):

  cdef list new_polys = [ R(0) ] * M.nrows()
  cdef int i, j, k

  for i from 0 <= i < M.nrows():
    for j from 0 <= j < M.rows[i].num_nonzero:
      k = M.rows[i].positions[j]
      new_polys[i] += M.rows[i].entries[j] * monomial_list[k]

  cdef MPolynomial_libsingular f, g

  cdef set previous_lms = set([ g.lm() for g in reducers ])

  for f in new_polys:
    if f != 0 and f.lm() not in previous_lms:
      G.append(clothed_polynomial(f, 1))

@cython.profile(True)
cpdef reduce_F4 (list L, set todo, list G):
  '''
  Build a matrix from a list of pairs L and reduces w.r.t. G.
  For now, I am assuming the normal selection strategy - this is relevant
  here because it means I don't have to worry about updating sugars.
  Faugère reports that the normal strategy works better in his original F4 paper.
  '''
  init_time = time.time()
  cdef Matrix_modn_sparse M, Mred
  cdef list reducers
  cdef list monomials
  cdef MPolynomialRing_libsingular R = L[0].parent()

  M, reducers, monomials = symbolic_preprocessing(L, todo, G) #Build the matrix
  statistics.inc_preprocessing_time(time.time() - init_time)

  before_red = time.time()
  Mred = M.rref() #Do row reduction
  statistics.inc_matrix_time(time.time() - before_red)

  before_add = time.time()
  add_polynomials_from_matrix(Mred, reducers, monomials, R, G)
  statistics.inc_addpolys_time(time.time() - before_add)

  statistics.inc_reduction_time(time.time() - init_time)

@cython.profile(True)
cpdef tuple select_pairs_normal_F4 (list P):
  '''
  Select pairs to be reduced in F4 using the normal strategy, that is,
  picks the critical pairs with minimal degree.
  '''

  cdef list L = []
  cdef set todo = set()

  if not P: #P is already empty
    return L, todo

  cdef int min_deg = P[0][2] #Degree of the first element of the S-polynomial
  #queue. Note that this works because we assume P is already sorted.

  cdef int i = 0
  cdef MPolynomial_libsingular f, g, tf, tg, tfg, s1, s2
  cdef MPolynomialRing_libsingular R
  while P and P[0][2] == min_deg:

    #compute both "branches" of the S-polynomial and add to the reducer list
    f = P[0][0].value()
    g = P[0][1].value()
    tfg = P[0][len(P[0]) - 1]

    if g == 0:

      s1 = f
      L.append(s1)
      todo.update(s1.monomials()[1:])

    else:

      tf = f.lm()
      tg = g.lm()
      R = tf.parent()
      s1 = R.monomial_quotient(tfg, tf) * f
      s2 = R.monomial_quotient(tfg, tg) * g

      #taking reducers from these polynomials
      #slicing off the first monomial of their lists (i.e., their lms)
      todo.update(s1.monomials()[1:], s2.monomials()[1:])

      L.append(s1)
      L.append(s2)

    P.pop(0)

  return L, todo

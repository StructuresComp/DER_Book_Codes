/**
 * \file PetscMatrix.inl
 *
 * \author miklos@cs.columbia.edu
 * \date 09/08/2009
 */

inline PetscMatrix::PetscMatrix(int s)
  : MatrixBase(s, s)
  , m_v(NULL)
  , m_w(NULL)
{
  setup(1);
}

inline PetscMatrix::PetscMatrix(int r, int c, int nnz)
  : MatrixBase(r, c)
  , m_v(NULL)
  , m_w(NULL)
{
  setup(nnz);
}

inline PetscMatrix::PetscMatrix(const PetscMatrix& M)
  : MatrixBase(0, 0)
  , m_v(NULL)
  , m_w(NULL)
{
  PetscInt m, n;
  MatGetSize(M.m_M, &m, &n);
  m_rows = m;
  m_cols = n;
  MatDuplicate(M.m_M, MAT_COPY_VALUES, &m_M);
  MatSetOption(m_M, MAT_USE_INODES, PETSC_FALSE);
  MatSetFromOptions(m_M);
  PetscUtils::createPetscVector(m_v, rows());
  PetscUtils::createPetscVector(m_w, cols());
}

inline PetscMatrix::~PetscMatrix()
{
  MatDestroy(m_M);
  if (m_v != NULL) VecDestroy(m_v);
  if (m_w != NULL) VecDestroy(m_w);
}

inline void PetscMatrix::setup(int nnz)
{
  MatCreateSeqAIJ(PETSC_COMM_SELF, m_rows, m_cols, nnz, PETSC_NULL, &m_M);
  MatSetOption(m_M, MAT_USE_INODES, PETSC_FALSE);
  MatSetFromOptions(m_M);
  PetscUtils::createPetscVector(m_v, rows());
  PetscUtils::createPetscVector(m_w, cols());
}

inline Scalar PetscMatrix::operator() (int r, int c) const
{
  PetscInt m = r;
  PetscInt n = c;
  PetscScalar val = 0;
  int ierr = MatGetValues(m_M, 1, &m, 1, &n, &val);
  return (Scalar) val;
}

inline int PetscMatrix::set(int r, int c, Scalar val)
{
  int ierr = MatSetValue(m_M, r, c, val, INSERT_VALUES);
  CHKERRQ(ierr);
  return 0;
}

inline int PetscMatrix::add(int r, int c, Scalar val)
{
  int ierr = MatSetValue(m_M, r, c, val, ADD_VALUES);
  CHKERRQ(ierr);
  return 0;
}

inline int PetscMatrix::add(const IntArray& rowIdx, const IntArray& colIdx,
                            const MatXd& values)
{
  int ierr;
  ierr = MatSetValues(m_M, rowIdx.size(), &(rowIdx[0]),
                      colIdx.size(), &(colIdx[0]),
                      values.data(), ADD_VALUES);
  CHKERRQ(ierr);
  return 0;
}

inline int PetscMatrix::add(const IndexArray& rowIdx, const IndexArray& colIdx,
                            const MatXd& values)
{
  int ierr;
  ierr = MatSetValues(m_M, rowIdx.size(), rowIdx.data(),
                      colIdx.size(), colIdx.data(),
                      values.data(), ADD_VALUES);
  CHKERRQ(ierr);
  return 0;
}

inline int PetscMatrix::scale(Scalar val)
{
  int ierr = MatScale(m_M, val);
  CHKERRQ(ierr);
  return 0;
}

inline int PetscMatrix::multiply(VecXd& y, Scalar s, const VecXd& x) const
{
  assert(cols() == x.size());
  assert(rows() == y.size());

  PetscUtils::copyToPetscVector(m_v, x);

  MatMult(m_M, m_v, m_w);

  PetscUtils::addFromPetscVector(y, s, m_w);

  return 0;
}
/*

  int PetscMatrix::diagonalScale(VecXd& a)
  {
  assert(m_rows == a.size());

  PetscUtils::copyPetscVec(a, m_w);
  int ierr = MatDiagonalScale(m_M, m_w, PETSC_NULL); CHKERRQ(ierr);

  return 0;
  }



  int PetscMatrix::shift(Scalar a)
  {
  int ierr = MatShift(m_M, a); CHKERRQ(ierr);
  return 0;
  }



  int PetscMatrix::multiply(const PetscMatrix& A, const PetscMatrix& B, PetscMatrix*& C)
  {
  int ierr;

  if (C == NULL) {
  C = new PetscMatrix(A.m_rows, B.m_cols);
  ierr = MatDestroy(C->m_M); CHKERRQ(ierr);
  ierr = MatMatMultSymbolic(A.m_M, B.m_M, 1, &(C->m_M)); CHKERRQ(ierr);
  }

  ierr = MatMatMultNumeric(A.m_M, B.m_M, C->m_M); CHKERRQ(ierr);

  return 0;
  }



  int PetscMatrix::PtA(const PetscMatrix& P, const PetscMatrix& A, PetscMatrix*& C)
  {
  int ierr;

  if (C == NULL) {
  C = new PetscMatrix(P.m_cols, A.m_cols);
  ierr = MatDestroy(C->m_M); CHKERRQ(ierr);
  ierr = MatMatMultTranspose(P.m_M, A.m_M, MAT_INITIALm_MATRIX, 1, &(C->m_M));
  CHKERRQ(ierr);
  }

  ierr = MatMatMultTranspose(P.m_M, A.m_M, MAT_REUSEm_MATRIX, 1, &(C->m_M));
  CHKERRQ(ierr);

  return 0;
  }



  int PetscMatrix::Pta(const PetscMatrix& P, const VecXd& a, VecXd& b)
  {
  assert(0);
  assert(P.m_rows == a.size());
  assert(P.m_cols == b.size());

  // set up Petsc vectors
  Vec v; PetscUtils::copyPetscVec(a, v);
  Vec w; PetscUtils::createPetscVec(b, w);

  // multiply
  MatMultTranspose(P.m_M, v, w);

  // add the result to w
  Scalar* wDat;  int ierr = VecGetArray(w, &wDat);  CHKERRQ(ierr);
  for (int i = 0; i < b.size(); i++) b[i] += wDat[i];
  ierr = VecRestoreArray(w, &wDat);  CHKERRQ(ierr);

  // clean up memory
  VecDestroy(v); VecDestroy(w);

  return 0;
  }
*/
inline int PetscMatrix::setZero()
{
  int ierr = MatZeroEntries(m_M);
  CHKERRQ(ierr);
  return 0;
}

inline int PetscMatrix::zeroRows(const IntArray& idx, Scalar diag)
{
  MatSetOption(m_M, MAT_KEEP_ZEROED_ROWS, PETSC_TRUE);
  int ierr = MatZeroRows(m_M, idx.size(), &idx[0], diag);
  CHKERRQ(ierr);
  return 0;
}
/*
  int PetscMatrix::assemble()
  {
  int ierr = MatAssemblyBegin(m_M, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_M, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  return 0;
  }
*/

inline int PetscMatrix::finalize()
{
  int ierr = MatAssemblyBegin(m_M, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m_M, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatSetOption(m_M, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);
  CHKERRQ(ierr);
  ierr = MatSetOption(m_M, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  CHKERRQ(ierr);
  return 0;
}

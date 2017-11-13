/**
 * \file PetscMatrix.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/08/2009
 */

#ifndef PETSCMATRIX_HH
#define PETSCMATRIX_HH

#include <petscmat.h>
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Petsc/PetscUtils.hh"

namespace BASim {

class PetscMatrix : public MatrixBase
{
public:

  PetscMatrix(int s);
  PetscMatrix(int r, int c, int nnz = 1);
  PetscMatrix(const PetscMatrix& M);
  ~PetscMatrix();

  Scalar operator() (int i, int j) const;
  int set(int r, int c, Scalar val);
  int add(int r, int c, Scalar val);
  int add(const IntArray& rowIdx, const IntArray& colIdx, const MatXd& values);
  int add(const IndexArray& rowIdx, const IndexArray& colIdx,
          const MatXd& values);
  int scale(Scalar val);
  int setZero();
  int zeroRows(const IntArray& idx, Scalar diag = 1.0);
  int multiply(VecXd& y, Scalar s, const VecXd& x) const;

  /*
  template <class IndexArray, class ValueMatrix>
  int addValues(IndexArray& rowIdx, IndexArray& colIdx, ValueMatrix& vals);

  // b = r*M*a
  int multiply(Scalar r, const VecXd& a, VecXd& b);

  int diagonalScale(VecXd &a);

  int shift(Scalar a);

  // C = A * B
  static int multiply(const PetscMatrix& A, const PetscMatrix& B, PetscMatrix*& C);

  // C = P^T * A
  static int PtA(const PetscMatrix& P, const PetscMatrix& A, PetscMatrix*& C);

  // b = P^T * a
  static int Pta(const PetscMatrix& P, const VecXd& a, VecXd& b);

  Scalar get(int r, int c) const;
  */

  //int assemble();

  int finalize();

  const Mat& getPetscMatrix() const { return m_M; }
  Mat& getPetscMatrix() { return m_M; }

protected:

  void setup(int nnz);

  Mat m_M;       // the underlying Petsc matrix
  Vec m_v, m_w;  // m_v is size rows(), m_w is size cols()
};

#include "PetscMatrix.inl"

} // namespace BASim

#endif // PETSCMATRIX_HH

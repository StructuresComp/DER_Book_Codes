/**
 * \file MKLLinearSolver.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 11/29/2009
 */

#include "MKLLinearSolver.hh"
#include "BASim/src/Core/Util.hh"
#include "BASim/src/Math/BandMatrix.hh"
#include "LAPACKPrototypes.h"

using namespace std;

namespace BASim {

MKLLinearSolver::MKLLinearSolver(BandMatrix& A)
  : LinearSolverBase(A)
{
  assert(A.rows() == A.cols());
  int n = A.rows();
  int kl = smart_cast<BandMatrix&>(A).kl();
  int ku = smart_cast<BandMatrix&>(A).ku();
  ipiv = new int[n];

  // ab holds the entries of the matrix. Space must be made for an
  // additional kl super-diagonals for LU factorization
  ab = new double[(2 * kl + ku + 1) * n];
}

MKLLinearSolver::~MKLLinearSolver()
{
  delete [] ipiv;
  delete [] ab;
}

static inline void convert(double* ab, BandMatrix& A, int kl, int ku, int n)
{
  int NUMROWS = 2 * kl + ku + 1;
  for (int j = 0; j < n; ++j) {
    for (int i = std::max(0, j - ku); i < std::min(n, j + kl + 1); ++i) {
      int row = kl + ku + i - j;
      int col = j;
      int offset = row + col * NUMROWS;
      ab[offset] = A(i, j);
    }
  }
}

int MKLLinearSolver::solve(VecXd& x, const VecXd& b)
{
  int n, kl, ku, nrhs, ldab, ldb, info;

  BandMatrix& A = smart_cast<BandMatrix&>(m_A);
  n = A.rows();
  kl = A.kl();
  ku = A.ku();
  nrhs = 1;
  ldab = 2 * kl + ku + 1;
  ldb = n;
  convert(ab, A, kl, ku, n);
  x = b;

  dgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, x.data(), &ldb, &info);

  // check return value for errors
  if (info < 0) {
    std::cerr << "Error in parameter " << -info << " in call to dgbsv"
              << std::endl;
    return -1;
  } else if (info > 0) {
    std::cerr << "Factor U is singular (" << info << std::endl;
    return -1;
  }
  
  assert( info == 0 );
  
  // Check the inf norm of the residual
  #ifdef DEBUG
    VecXd residual(x.size());
    residual.setZero();
    A.multiply(residual,1.0,x);
    residual -= b;
    double infnorm = fabs(residual.maxCoeff());
    if( infnorm > 1.0e-6 )
    {
      std::cout << "\033[31;1mWARNING IN MKLLinearSolver:\033[m Large residual detected. ||residual||_{inf} = " << infnorm << std::endl;
      return -1;
    }
  #endif
  
  return 0;
}

} // namespace BASim

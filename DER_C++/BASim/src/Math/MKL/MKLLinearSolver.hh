/**
 * \file MKLLinearSolver.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 11/29/2009
 */

#ifndef MKLLINEARSOLVER_HH
#define MKLLINEARSOLVER_HH

#include "BASim/src/Math/LinearSolverBase.hh"

namespace BASim {

class BandMatrix;

/** LU-based linear solver for band matrices. */
class MKLLinearSolver : public LinearSolverBase
{
public:

  MKLLinearSolver(BandMatrix& A);
  virtual ~MKLLinearSolver();

  virtual int solve(VecXd& x, const VecXd& b);

private:

  int* ipiv;
  double* ab;
};

} // namespace BASim

#endif // MKLLINEARSOLVER_HH

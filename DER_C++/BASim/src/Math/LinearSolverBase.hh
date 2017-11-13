/**
 * \file LinearSolverBase.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/08/2009
 */

#ifndef LINEARSOLVERBASE_HH
#define LINEARSOLVERBASE_HH

#include "BASim/src/Math/MatrixBase.hh"

namespace BASim {

/** Interface for linear solvers. */
class LinearSolverBase
{
public:

  LinearSolverBase(MatrixBase& A)
    : m_A(A)
  {}

  virtual ~LinearSolverBase() {}

  /**
   * Solves the equation \f$Ax=b\f$ for \f$x\f$, given \f$A\f$ and
   * \f$b\f$.
  */
  virtual int solve(VecXd& x, const VecXd& b) = 0;

protected:

  MatrixBase& m_A;
};

} // namespace BASim

#endif // LINEARSOLVERBASE_HH

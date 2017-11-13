/**
 * \file DiagonalPreconditioner.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 11/16/2009
 */

#ifndef DIAGONALPRECONDITIONER_HH
#define DIAGONALPRECONDITIONER_HH

#include "BASim/src/Math/Preconditioner.hh"
#include "BASim/src/Math/MatrixBase.hh"

namespace BASim {

class DiagonalPreconditioner : public Preconditioner
{
public:

  DiagonalPreconditioner(const MatrixBase& M)
    : diagonals(M.rows())
  {
    for (int i = 0; i < M.rows(); ++i) {
      Scalar val = M(i, i);
      if (val == 0) diagonals(i) = 1;
      else diagonals(i) = 1.0 / (val > 0 ? val : -val);
    }
  }

  virtual void apply(VecXd& w, const VecXd& v)
  {
    assert(w.size() == diagonals.size());
    assert(v.size() == diagonals.size());

    int s = diagonals.size();
    for (int i = 0; i < s; ++i) { w(i) = diagonals(i) * v(i); }
  }

protected:

  VecXd diagonals;
};

} // namespace BASim

#endif // DIAGONALPRECONDITIONER_HH

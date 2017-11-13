/**
 * \file Preconditioner.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 11/16/2009
 */

#ifndef PRECONDITIONER_HH
#define PRECONDITIONER_HH

#include "BASim/src/Core/Definitions.hh"

namespace BASim {

class Preconditioner
{
public:
  virtual ~Preconditioner() {}
  virtual void apply(VecXd& w, const VecXd& v) = 0;
};

} // namespace BASim

#endif // PRECONDITIONER_HH

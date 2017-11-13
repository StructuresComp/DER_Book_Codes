/**
 * \file Force.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/14/2009
 */

#ifndef FORCE_HH
#define FORCE_HH

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Math/MatrixBase.hh"

namespace BASim {

/** */
class Force
{
public:

  Force() {}

  virtual ~Force() {}

  virtual void globalForce(VecXd& force) = 0;
  virtual void globalJacobian(MatrixBase& Jacobian) = 0;
};

} // namespace BASim

#endif // FORCE_HH

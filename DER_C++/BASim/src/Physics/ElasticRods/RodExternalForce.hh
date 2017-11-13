/**
 * \file RodExternalForce.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/04/2009
 */

#ifndef RODEXTERNALFORCE_HH
#define RODEXTERNALFORCE_HH

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"

namespace BASim {

/** Base class for external forces applied to a rod. */
class RodExternalForce
{
public:

  RodExternalForce(bool implicit = true) : m_implicit(implicit) {}

  virtual ~RodExternalForce() {}

//  virtual void computeForce(const ElasticRod& rod, VecXd& force) = 0;
  virtual void computeForce(ElasticRod& rod, VecXd& force) = 0;
  virtual void computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) = 0;
  virtual void computeForceDV(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) = 0;

  bool isImplicit() const { return m_implicit; }
  void setImplicit(bool implicit) { m_implicit = implicit; }

protected:

  bool m_implicit;
};

} // namespace BASim

#endif // RODEXTERNALFORCE_HH

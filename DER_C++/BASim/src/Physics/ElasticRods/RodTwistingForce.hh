/**
 * \file RodTwistingForce.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/06/2009
 */

#ifndef RODTWISTINGFORCE_HH
#define RODTWISTINGFORCE_HH

#include "BASim/src/Physics/ElasticRods/RodForce.hh"
#include "BASim/src/Physics/ElasticRods/VertexStencil.hh"

namespace BASim {

class RodTwistingForce : public RodForceT<VertexStencil>
{
public:

  typedef Eigen::Matrix<Scalar, 11, 1> ElementForce;
  typedef Eigen::Matrix<Scalar, 11, 11> ElementJacobian;

  RodTwistingForce(ElasticRod& rod);

  virtual Scalar globalEnergy();
  virtual void globalForce(VecXd& force);
  virtual void globalJacobian(int baseidx, Scalar scale, MatrixBase& J);

  Scalar localEnergy(const vertex_handle& vh);
  void localForce(ElementForce& force, const vertex_handle& vh);
  void localJacobian(ElementJacobian& Jacobian, const vertex_handle& vh);

  const Scalar& getKt(const vertex_handle& vh) const;
  void setKt(const vertex_handle& vh, const Scalar& kt);

  const Scalar& getTwist(const vertex_handle& vh) const;
  void setTwist(const vertex_handle& vh, const Scalar& twist);

  const Scalar& getUndeformedTwist(const vertex_handle& vh) const;
  void setUndeformedTwist(const vertex_handle& vh,
                          const Scalar& undeformedTwist);

  const Scalar& getRefVertexLength(const vertex_handle& vh) const;
  void setRefVertexLength(const vertex_handle& vh, const Scalar& length);

  virtual void updateProperties();
  virtual void updateStiffness();
  virtual void updateUndeformedStrain();
  virtual void updateReferenceDomain();

protected:

#ifdef TEST_ROD_TWISTING
  void testEnergy(const Scalar& energy, const vertex_handle& vh) const;
  void testForce(const ElementForce& force, const vertex_handle& vh) const;
  void testJacobian(const ElementJacobian& Jacobian,
                    const vertex_handle& vh) const;
#endif // TEST_ROD_TWISTING

  VPropHandle<Scalar> m_kt;              ///< twist stiffness
  VPropHandle<Scalar> m_twist;           ///< twist at a vertex
  VPropHandle<Scalar> m_undeformedTwist; ///< undeformed twist
  VPropHandle<Scalar> m_refVertexLength; ///< length of domain of integration
};

} // namespace BASim

#endif // RODTWISTINGFORCE_HH

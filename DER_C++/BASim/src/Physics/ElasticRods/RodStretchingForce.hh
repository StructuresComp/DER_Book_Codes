/**
 * \file RodStretchingForce.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/01/2009
 */

#ifndef RODSTRETCHINGFORCE_HH
#define RODSTRETCHINGFORCE_HH

#include "BASim/src/Physics/ElasticRods/RodForce.hh"
#include "BASim/src/Physics/ElasticRods/EdgeStencil.hh"

namespace BASim {

class SpringDofStruct
{
public:
  Vec3d x[2];           ///< position of two endpoints
  Vec3d edge;           ///< vector from x[0] to x[1]
  Vec3d tangent;        ///< normalized vector from x[0] to x[1]
  Scalar currLength;    ///< current length
  Scalar restLength;    ///< rest length
  Scalar stiffness;     ///< spring stiffness
};

/** This class implements the stretching force for an elastic rod. */
class RodStretchingForce : public RodForceT<EdgeStencil>
{
public:

  typedef Eigen::Matrix<Scalar, 6, 1> ElementForce;
  typedef Eigen::Matrix<Scalar, 6, 6> ElementJacobian;

  RodStretchingForce(ElasticRod& rod);

  void gatherDofs(SpringDofStruct& dofs, const edge_handle& eh);

  virtual Scalar globalEnergy();
  virtual void globalForce(VecXd& force);
  virtual void globalJacobian(int baseidx, Scalar scale, MatrixBase& Jacobian);

  Scalar elementEnergy(const edge_handle& eh);
  void elementForce(ElementForce& force, const SpringDofStruct& dofs);
  void elementForce(ElementForce& force, const edge_handle& eh);
  void elementJacobian(ElementJacobian& Jacobian, const edge_handle& eh);

  const Scalar& getKs(const edge_handle& eh) const;
  void setKs(const edge_handle& eh, const Scalar& ks);

  const Scalar& getRefLength(const edge_handle& eh) const;
  void setRefLength(const edge_handle& eh, const Scalar& length);

  void updateStiffness();
  void updateUndeformedStrain();
  
	virtual void updatePropertiesForNewVertex();  

protected:

#ifdef TEST_ROD_STRETCHING
  void testEnergy(const Scalar& energy, const edge_handle& eh) const;
  void testForce(const ElementForce& force, const edge_handle& eh) const;
  void testJacobian(const ElementJacobian& Jacobian,
                    const edge_handle& eh) const;
#endif // TEST_ROD_STRETCHING

  EPropHandle<Scalar> m_ks;
  EPropHandle<Scalar> m_refLength;
};

} // namespace BASim

#endif // RODSTRETCHINGFORCE_HH

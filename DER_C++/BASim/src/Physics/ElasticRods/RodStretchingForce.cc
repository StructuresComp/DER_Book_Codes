/**
 * \file RodStretchingForce.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/01/2009
 */

#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
#include "BASim/src/Math/Math.hh"

#ifdef TEST_ROD_STRETCHING
#include "BASim/src/Physics/ElasticRods/Tests/RodStretchingTest.hh"
#endif // TEST_ROD_STRETCHING

using namespace std;

namespace BASim {

RodStretchingForce::RodStretchingForce(ElasticRod& rod)
  : RodForceT<EdgeStencil>(rod, "RodStretchingForce")
{
  m_rod.add_property(m_ks, "stretching stiffness");
  m_rod.add_property(m_refLength, "stretching ref length");

  updateStiffness();
  updateUndeformedStrain();
}

void RodStretchingForce::updatePropertiesForNewVertex() {
  iterator end = m_stencil.end();
  --end;
  
  m_stencil = end;

	{
    edge_handle& eh = m_stencil.handle();
    Scalar a = m_rod.radiusA(eh);
    Scalar b = m_rod.radiusB(eh);
    Scalar E = m_rod.getYoungsModulus();
    if (viscous()) {
      E = 3 * m_rod.getViscosity() / m_rod.getTimeStep();
    }
    Scalar A = a * b;    
    if (m_rod.getIsHollow()) A = A - m_rod.getInnerRadius() * m_rod.getInnerRadius();    
    setKs(eh, E * M_PI * A);
  }

	{  
    edge_handle& eh = m_stencil.handle();
    setRefLength(eh, m_rod.getEdgeLength(eh));
  }
  

}

void RodStretchingForce::gatherDofs(SpringDofStruct& dofs,
                                    const edge_handle& eh)
{
  dofs.x[0] = m_rod.getFromVertex(eh);
  dofs.x[1] = m_rod.getToVertex(eh);
  dofs.edge = m_rod.getEdge(eh);
  dofs.tangent = m_rod.getTangent(eh);
  dofs.currLength = m_rod.getEdgeLength(eh);
  dofs.restLength = getRefLength(eh);
  dofs.stiffness = getKs(eh);
}

Scalar RodStretchingForce::globalEnergy()
{
  Scalar energy = 0;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    edge_handle& eh = m_stencil.handle();
    Scalar localEnergy = elementEnergy(eh);
    energy += localEnergy;

#ifdef TEST_ROD_STRETCHING
    testEnergy(localEnergy, eh);
#endif // TEST_ROD_STRETCHING

  }
  return energy;
}

Scalar RodStretchingForce::elementEnergy(const edge_handle& eh)
{
  Scalar ks = getKs(eh);
  if (ks == 0.0) return 0;

  Scalar refLength = getRefLength(eh);
  Scalar len = m_rod.getEdgeLength(eh);

  return ks / 2.0 * square(len / refLength - 1.0) * refLength;
}

void RodStretchingForce::globalForce(VecXd& force)
{
  IndexArray indices;
  ElementForce localForce;
  SpringDofStruct dofs;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    edge_handle& eh = m_stencil.handle();
    //elementForce(localForce, eh);
    gatherDofs(dofs, eh);
    elementForce(localForce, dofs);
    m_stencil.indices(indices);
    for (int i = 0; i < indices.size(); ++i)
      force(indices(i)) += localForce(i);

#ifdef TEST_ROD_STRETCHING
    testForce(localForce, eh);
#endif // TEST_ROD_STRETCHING
  }

#ifdef TEST_ROD_STRETCHING
  globalEnergy();
#endif // TEST_ROD_STRETCHING
}

void RodStretchingForce::elementForce(ElementForce& force,
                                      const SpringDofStruct& dofs)
{
  Vec3d f =
    dofs.stiffness * (dofs.currLength / dofs.restLength - 1.0) * dofs.tangent;
  force.segment(0, 3) =  f;
  force.segment(3, 3) = -f;
}

void RodStretchingForce::elementForce(ElementForce& force,
                                      const edge_handle& eh)
{
  Scalar ks = getKs(eh);
  if (ks == 0.0) return;

  Scalar refLength = getRefLength(eh);
  Scalar len = m_rod.getEdgeLength(eh);
  const Vec3d& tangent = m_rod.getTangent(eh);

  Vec3d f = ks * (len / refLength - 1.0) * tangent;
  force.segment(0, 3) =  f;
  force.segment(3, 3) = -f;
}

void RodStretchingForce::globalJacobian(int baseidx, Scalar scale, MatrixBase& Jacobian)
{
  IndexArray indices;
  ElementJacobian localJ;
  MatXd adder;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    edge_handle& eh = m_stencil.handle();
    localJ.setZero();
    elementJacobian(localJ, eh);
    adder = localJ;
    adder *= scale;
    m_stencil.indices(indices);
    for( int i = 0; i < (int) indices.size(); ++i ) indices(i) += baseidx;
    Jacobian.add(indices, indices, adder);

#ifdef TEST_ROD_STRETCHING
    testJacobian(localJ, eh);
#endif // TEST_ROD_STRETCHING
  }
}

void RodStretchingForce::elementJacobian(ElementJacobian& Jacobian,
                                         const edge_handle& eh)
{
  Scalar ks = getKs(eh);
  if (ks == 0.0) return;

  const Vec3d& e = m_rod.getEdge(eh);
  Scalar len = m_rod.getEdgeLength(eh);
  Scalar refLength = getRefLength(eh);
  Mat3d M = ks * ( (1.0 / refLength - 1.0 / len) * Mat3d::Identity()
                   + 1.0 / len * outerProd(e,e) / square(len) );

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      Jacobian(    i, j) = Jacobian(3 + i, 3 + j) = -M(i, j);
      Jacobian(3 + i, j) = Jacobian(    i, 3 + j) =  M(i, j);
    }
  }
}

const Scalar& RodStretchingForce::getKs(const edge_handle& eh) const
{
  return m_rod.property(m_ks)[eh];
}

void RodStretchingForce::setKs(const edge_handle& eh, const Scalar& ks)
{
  m_rod.property(m_ks)[eh] = ks;
}

const Scalar& RodStretchingForce::getRefLength(const edge_handle& eh) const
{
  return m_rod.property(m_refLength)[eh];
}

void RodStretchingForce::setRefLength(const edge_handle& eh,
                                      const Scalar& length)
{
  m_rod.property(m_refLength)[eh] = length;
}

void RodStretchingForce::updateStiffness()
{
  Scalar E = m_rod.getYoungsModulus();
  if (viscous()) {
    E = 3 * m_rod.getViscosity() / m_rod.getTimeStep();
  }

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    edge_handle& eh = m_stencil.handle();
    Scalar a = m_rod.radiusA(eh);
    Scalar b = m_rod.radiusB(eh);
    Scalar A = a * b;    
    if (m_rod.getIsHollow()) A = A - m_rod.getInnerRadius() * m_rod.getInnerRadius();    
    setKs(eh, E * M_PI * A);
  }
}

void RodStretchingForce::updateUndeformedStrain()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    edge_handle& eh = m_stencil.handle();
    setRefLength(eh, m_rod.getEdgeLength(eh));
  }
}

#ifdef TEST_ROD_STRETCHING

void RodStretchingForce::testEnergy(const Scalar& energy,
                                    const edge_handle& eh) const
{
  Scalar ks = getKs(eh);
  Scalar refLength = getRefLength(eh);
  Scalar len = m_rod.getEdgeLength(eh);
  const Vec3d& x0 = m_rod.getFromVertex(eh);
  const Vec3d& x1 = m_rod.getToVertex(eh);
  Scalar mathEnergy;
  rodStretchingEnergyTest(mathEnergy, energy, x0, x1, ks, refLength);
}

void RodStretchingForce::testForce(const ElementForce& force,
                                   const edge_handle& eh) const
{
  Scalar ks = getKs(eh);
  Scalar refLength = getRefLength(eh);
  Scalar len = m_rod.getEdgeLength(eh);
  const Vec3d& tangent = m_rod.getTangent(eh);
  const Vec3d& x0 = m_rod.getFromVertex(eh);
  const Vec3d& x1 = m_rod.getToVertex(eh);
  ElementForce mathForce;
  rodStretchingForceTest(mathForce, force, x0, x1, ks, refLength);
}

void RodStretchingForce::testJacobian(const ElementJacobian& Jacobian,
                                      const edge_handle& eh) const
{
  const Vec3d& x0 = m_rod.getFromVertex(eh);
  const Vec3d& x1 = m_rod.getToVertex(eh);
  Scalar ks = getKs(eh);
  Scalar refLength = getRefLength(eh);
  ElementJacobian mathJ;
  rodStretchingJacobianTest(mathJ, Jacobian, x0, x1, ks, refLength);
}

#endif // TEST_ROD_STRETCHING

} // namespace BASim

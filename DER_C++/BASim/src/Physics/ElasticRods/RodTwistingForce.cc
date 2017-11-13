/**
 * \file RodTwistingForce.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/06/2009
 */

#include "BASim/src/Physics/ElasticRods/RodTwistingForce.hh"
#include "BASim/src/Math/Math.hh"

#ifdef TEST_ROD_TWISTING
#include "BASim/src/Physics/ElasticRods/Tests/RodTwistingTest.hh"
#endif // TEST_ROD_TWISTING

namespace BASim {

RodTwistingForce::RodTwistingForce(ElasticRod& rod)
  : RodForceT<VertexStencil>(rod, "RodTwistingForce")
{
  m_rod.add_property(m_kt, "twist stiffness");
  m_rod.add_property(m_twist, "twist");
  m_rod.add_property(m_undeformedTwist, "undeformed twist");
  m_rod.add_property(m_refVertexLength, "twisting ref length");

  updateProperties();
  updateUndeformedStrain();
  updateStiffness();
  updateReferenceDomain();
}

void RodTwistingForce::updateProperties()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    edge_handle eh0 = m_stencil.inEdge();
    edge_handle eh1 = m_stencil.outEdge();
    setTwist(vh, m_rod.getReferenceTwist(vh)
             + m_rod.getTheta(eh1) - m_rod.getTheta(eh0));
  }
}

void RodTwistingForce::updateStiffness()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    edge_handle eh0 = m_stencil.inEdge();
    edge_handle eh1 = m_stencil.outEdge();
    Scalar G = m_rod.getShearModulus();
    Scalar a = (m_rod.radiusA(eh0) + m_rod.radiusA(eh1)) / 2.0;
    Scalar b = (m_rod.radiusB(eh0) + m_rod.radiusB(eh1)) / 2.0;
    setKt(vh, M_PI * (pow(a, 4) + pow(b, 4)) * G / 4.0);
  }
}

void RodTwistingForce::updateReferenceDomain()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    edge_handle eh0 = m_stencil.inEdge();
    edge_handle eh1 = m_stencil.outEdge();
    Scalar len = (m_rod.getEdgeLength(eh0) + m_rod.getEdgeLength(eh1)) / 2.0;
    setRefVertexLength(vh, len);
  }
}

Scalar RodTwistingForce::globalEnergy()
{
  Scalar energy = 0;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    Scalar localE = localEnergy(vh);
    energy += localE;

#ifdef TEST_ROD_TWISTING
    testEnergy(localE, vh);
#endif // TEST_ROD_TWISTING
  }

  return energy;
}

Scalar RodTwistingForce::localEnergy(const vertex_handle& vh)
{
  Scalar kt = getKt(vh);
  Scalar undefTwist = getUndeformedTwist(vh);
  Scalar len = getRefVertexLength(vh);
  return kt / (2.0 * len) * square(getTwist(vh) - undefTwist);
}

void RodTwistingForce::globalForce(VecXd& force)
{
  ElementForce f;
  IndexArray indices;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    localForce(f, vh);
    m_stencil.indices(indices);
    for (int j = 0; j < f.size(); ++j) force(indices[j]) += f(j);

#ifdef TEST_ROD_TWISTING
    testForce(f, vh);
#endif // TEST_ROD_TWISTING
  }
#ifdef TEST_ROD_TWISTING
  globalEnergy();
#endif // TEST_ROD_TWISTING
}

void RodTwistingForce::localForce(ElementForce& force,
                                  const vertex_handle& vh)
{
  const Vec3d& kb = m_rod.getCurvatureBinormal(vh);
  Scalar twist = getTwist(vh);
  Scalar undefTwist = getUndeformedTwist(vh);
  Scalar kt = getKt(vh);
  Scalar len = getRefVertexLength(vh);
  Scalar len0 = m_rod.getEdgeLength(m_stencil.inEdge());
  Scalar len1 = m_rod.getEdgeLength(m_stencil.outEdge());

  force(3) =  kt/len * (twist - undefTwist);
  force(7) = -kt/len * (twist - undefTwist);

  for (int i = 0; i < 3; ++i) {
    force(i) = kt/len * (twist - undefTwist) / (2.0 * len0) * kb(i);
    force(8 + i) = -kt/len * (twist - undefTwist) / (2.0 * len1) * kb(i);
    force(4 + i) = -(force(i) + force(8 + i));
  }
}

void RodTwistingForce::globalJacobian(int baseidx, Scalar scale, MatrixBase& J)
{
  ElementJacobian localJ;
  MatXd adder;
  IndexArray indices;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    localJ.setZero();
    localJacobian(localJ, vh);
    m_stencil.indices(indices);
    for( int i = 0; i < (int) indices.size(); ++i ) indices[i] += baseidx;
    adder = localJ;
    adder *= scale;
    J.add(indices, indices, localJ);

#ifdef TEST_ROD_TWISTING
    testJacobian(localJ, vh);
#endif // TEST_ROD_TWISTING
  }
}

void RodTwistingForce::localJacobian(ElementJacobian& Jacobian,
                                     const vertex_handle& vh)
{
  const Vec3d& x0 = m_rod.getVertex(m_stencil.prevVert());
  const Vec3d& x1 = m_rod.getVertex(vh);
  const Vec3d& x2 = m_rod.getVertex(m_stencil.nextVert());

  const Vec3d& e0 = m_rod.getEdge(m_stencil.inEdge());
  const Vec3d& e1 = m_rod.getEdge(m_stencil.outEdge());
  const Vec3d& kb = m_rod.getCurvatureBinormal(vh);
  Scalar twist = getTwist(vh);
  Scalar len = getRefVertexLength(vh);
  Scalar kt = getKt(vh);
  Scalar undefTwist = getUndeformedTwist(vh);

  // d^2 E / (d\theta_i d\theta_j)
  Jacobian(3,3) = Jacobian(7,7) = -kt/len;
  Jacobian(3,7) = Jacobian(7,3) =  kt/len;

  // d^2 E / (d\x_i d\theta_j)
  Vec3d dEdx0 = -kt/(2 * len * e0.norm()) * kb;
  Vec3d dEdx2 =  kt/(2 * len * e1.norm()) * kb;
  Vec3d dEdx1 = -(dEdx0 + dEdx2);
  for (int i = 0; i < 3; ++i) {
    // x0,theta0
    Jacobian(i, 3) = Jacobian(3, i) =  dEdx0(i);
    // x0,theta1
    Jacobian(i, 7) = Jacobian(7, i) = -dEdx0(i);

    // x1,theta0
    Jacobian(4+i, 3) = Jacobian(3, 4+i) =  dEdx1(i);
    // x1,theta1
    Jacobian(4+i, 7) = Jacobian(7, 4+i) = -dEdx1(i);

    // x2,theta0
    Jacobian(8+i, 3) = Jacobian(3, 8+i) =  dEdx2(i);
    // x2,theta1
    Jacobian(8+i, 7) = Jacobian(7, 8+i) = -dEdx2(i);
  }

  // d^2 E / (d\x_i d\x_j)
  Mat3dArray Dkb = computeDkb(x0, x1, x2);

  Mat3d kbkb = outerProd(kb, kb);
  Mat3d kbe0 = outerProd(kb, e0);
  Mat3d kbe1 = outerProd(kb, e1);

  Mat3d term1 = -kt / (4 * len * e0.squaredNorm())
    * (kbkb - 2 * (twist - undefTwist) / e0.norm() * kbe0);
  Mat3d term2 = -kt / (4 * len * e0.norm() * e1.norm()) * kbkb;
  Mat3d term3 = -kt / (4 * len * e1.squaredNorm())
    * (kbkb - 2 * (twist - undefTwist) / e1.norm() * kbe1);
  Scalar t0 = -kt * (twist - undefTwist) / (2 * len * e0.norm());
  Scalar t1 = -kt * (twist - undefTwist) / (2 * len * e1.norm());
  Jacobian.block(0, 0, 3, 3) = term1 - t0 * Dkb[0];
  Jacobian.block(0, 4, 3, 3) = -term1 + term2 - t0 * Dkb[1];
  Jacobian.block(0, 8, 3, 3) = -term2 - t0 * Dkb[2];
  Jacobian.block(4, 0, 3, 3) = -term1 + term2 + (t0 - t1) * Dkb[0];
  Jacobian.block(4, 4, 3, 3) = term1 - 2.0 * term2 + term3 + (t0 - t1) * Dkb[1];
  Jacobian.block(4, 8, 3, 3) = term2 - term3 + (t0 - t1) * Dkb[2];
  Jacobian.block(8, 0, 3, 3) = -term2 + t1 * Dkb[0];
  Jacobian.block(8, 4, 3, 3) = term2 - term3 + t1 * Dkb[1];
  Jacobian.block(8, 8, 3, 3) = term3 + t1 * Dkb[2];
}

void RodTwistingForce::updateUndeformedStrain()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    setUndeformedTwist(vh, getTwist(vh));
  }
}

const Scalar& RodTwistingForce::getKt(const vertex_handle& vh) const
{
  return m_rod.property(m_kt)[vh];
}

void RodTwistingForce::setKt(const vertex_handle& vh, const Scalar& kt)
{
  m_rod.property(m_kt)[vh] = kt;
}

const Scalar& RodTwistingForce::getTwist(const vertex_handle& vh) const
{
  return m_rod.property(m_twist)[vh];
}

void RodTwistingForce::setTwist(const vertex_handle& vh, const Scalar& twist)
{
  m_rod.property(m_twist)[vh] = twist;
}

const Scalar&
RodTwistingForce::getUndeformedTwist(const vertex_handle& vh) const
{
  return m_rod.property(m_undeformedTwist)[vh];
}

void RodTwistingForce::setUndeformedTwist(const vertex_handle& vh,
                                          const Scalar& undeformedTwist)
{
  m_rod.property(m_undeformedTwist)[vh] = undeformedTwist;
}

const Scalar&
RodTwistingForce::getRefVertexLength(const vertex_handle& vh) const
{
  return m_rod.property(m_refVertexLength)[vh];
}

void RodTwistingForce::setRefVertexLength(const vertex_handle& vh,
                                          const Scalar& length)
{
  m_rod.property(m_refVertexLength)[vh] = length;
}

#ifdef TEST_ROD_TWISTING

void RodTwistingForce::testEnergy(const Scalar& energy,
                                  const vertex_handle& vh) const
{
  Scalar referenceTwist = m_rod.getReferenceTwist(vh);
  Scalar theta0 = m_rod.getTheta(m_stencil.inEdge());
  Scalar theta1 = m_rod.getTheta(m_stencil.outEdge());
  Scalar kt = getKt(vh);
  Scalar undefTwist = getUndeformedTwist(vh);
  Scalar len = getRefVertexLength(vh);
  Scalar mathEnergy;
  rodTwistingEnergyTest(mathEnergy, energy, referenceTwist, theta0, theta1,
                        undefTwist, kt, len);
}

void RodTwistingForce::testForce(const ElementForce& force,
                                 const vertex_handle& vh) const
{
  const Vec3d& x0 = m_rod.getVertex(m_stencil.prevVert());
  const Vec3d& x1 = m_rod.getVertex(vh);
  const Vec3d& x2 = m_rod.getVertex(m_stencil.nextVert());
  Scalar theta0 = m_rod.getTheta(m_stencil.inEdge());
  Scalar theta1 = m_rod.getTheta(m_stencil.outEdge());
  Scalar kt = getKt(vh);
  Scalar referenceTwist = m_rod.getReferenceTwist(vh);
  Scalar undefTwist = getUndeformedTwist(vh);
  Scalar len = getRefVertexLength(vh);
  ElementForce mathF;
  rodTwistingForceTest(mathF, force, x0, x1, x2, theta0, theta1,
                       referenceTwist, undefTwist, kt, len);
}

void RodTwistingForce::testJacobian(const ElementJacobian& Jacobian,
                                    const vertex_handle& vh) const
{
  const Vec3d& x0 = m_rod.getVertex(m_stencil.prevVert());
  const Vec3d& x1 = m_rod.getVertex(vh);
  const Vec3d& x2 = m_rod.getVertex(m_stencil.nextVert());
  Scalar theta0 = m_rod.getTheta(m_stencil.inEdge());
  Scalar theta1 = m_rod.getTheta(m_stencil.outEdge());
  Scalar kt = getKt(vh);
  Scalar referenceTwist = m_rod.getReferenceTwist(vh);
  Scalar undefTwist = getUndeformedTwist(vh);
  Scalar len = getRefVertexLength(vh);
  ElementJacobian mathJ;
  rodTwistingJacobianTest(mathJ, Jacobian, x0, x1, x2, theta0, theta1,
                          referenceTwist, undefTwist, kt, len);
}

#endif // TEST_ROD_TWISTING

} // namespace BASim

/**
 * \file RodBendingForce.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/06/2009
 */

#include "BASim/src/Physics/ElasticRods/RodBendingForce.hh"
#include "BASim/src/Math/Math.hh"

#ifdef TEST_ROD_BENDING
#include "BASim/src/Physics/ElasticRods/Tests/RodBendingTest.hh"
#endif // TEST_ROD_BENDING

namespace BASim {

RodBendingForce::RodBendingForce(ElasticRod& rod)
  : RodAnisoBending(rod, "RodBendingForce")
{
  m_rod.add_property(m_denominators, "denominators");
  m_rod.add_property(m_derivKb, "curvature binormal derivative");
  m_rod.add_property(m_edgeRotationMatrix, "rotation matrix from reference to material frame");
  updateProperties();
}

void RodBendingForce::updateProperties()
{
  RodAnisoBending::updateProperties();

  computeRotationMatrix();
  computeDenominators();
  computeDkb();

  firstDerivs = false;
  secondDerivs = false;
}

void RodBendingForce::computeRotationMatrix()
{
  for (int j = 0; j < m_rod.ne(); ++j) {
    const Scalar& theta = m_rod.getTheta(j);
    Mat2d& Rot = m_rod.property(m_edgeRotationMatrix)[j];
    Rot << -sin(theta), cos(theta), -cos(theta), -sin(theta);
  }
}

void RodBendingForce::computeDenominators()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    edge_handle eh0 = m_stencil.inEdge();
    edge_handle eh1 = m_stencil.outEdge();
    const Vec3d& e0 = m_rod.getEdge(eh0);
    const Vec3d& e1 = m_rod.getEdge(eh1);
    const Scalar& e0len = m_rod.getEdgeLength(eh0);
    const Scalar& e1len = m_rod.getEdgeLength(eh1);
    m_rod.property(m_denominators)[vh] = e0len * e1len + e0.dot(e1);
  }
}

void RodBendingForce::computeDxDenominator()
{
  m_DxDenominator.resize(m_rod.nv());

  for (int i = 1; i < m_rod.nv()-1; ++i) {
    const Vec3d& e0 = m_rod.getEdge(i - 1);
    const Vec3d& e1 = m_rod.getEdge(i);

    Vec9d& DxDenominator = m_DxDenominator[i];

    DxDenominator.segment(0,3) = -e1 - e1.norm() / e0.norm() * e0;
    DxDenominator.segment(6,3) =  e0 + e0.norm() / e1.norm() * e1;
    DxDenominator.segment(3,3) = -(DxDenominator.segment(0,3)+DxDenominator.segment(6,3));
  }
}

void RodBendingForce::computeDxDxDenominator()
{
  m_DxDxDenominator.resize(m_rod.nv());
  Mat3d Id = Mat3d::Identity();

  for (int k = 1; k < m_rod.nv()-1; ++k) {
    Mat9d& DxDxDenominator = m_DxDxDenominator[k];

    const Vec3d& e0 = m_rod.getEdge(k - 1);
    const Vec3d& e1 = m_rod.getEdge(k);

    Mat3d Id = Mat3d::Identity();
    Mat3d term1 = Id + outerProd(e1, e0) / (e0.norm() * e1.norm());
    Mat3d term2 = e0.norm() / e1.norm() * Id
      - e0.norm() / cube(e1.norm()) * outerProd(e1, e1);
    Mat3d term3 = Id + outerProd(e0, e1) / (e0.norm() * e1.norm());
    Mat3d term4 = e1.norm() / e0.norm() * Id
      - e1.norm() / cube(e0.norm()) * outerProd(e0, e0);

    DxDxDenominator.block(0, 0, 3, 3) = term4;
    DxDxDenominator.block(0, 3, 3, 3) = term3 - term4;
    DxDxDenominator.block(0, 6, 3, 3) = -term3;
    DxDxDenominator.block(3, 0, 3, 3) = term1 - term4;
    DxDxDenominator.block(3, 3, 3, 3) = -term1 + term2 - term3 + term4;
    DxDxDenominator.block(3, 6, 3, 3) = -term2 + term3;
    DxDxDenominator.block(6, 0, 3, 3) = -term1;
    DxDxDenominator.block(6, 3, 3, 3) = term1 - term2;
    DxDxDenominator.block(6, 6, 3, 3) = term2;
  }
}

void RodBendingForce::computeDkb()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    edge_handle eh0 = m_stencil.inEdge();
    edge_handle eh1 = m_stencil.outEdge();
    const Vec3d& t0 = m_rod.getTangent(eh0);
    const Vec3d& t1 = m_rod.getTangent(eh1);
    const Scalar& e0len = m_rod.getEdgeLength(eh0);
    const Scalar& e1len = m_rod.getEdgeLength(eh1);
    const Vec3d& kb = m_rod.getCurvatureBinormal(vh);

    Mat3dArray& Dkb = m_rod.property(m_derivKb)[vh];
    Dkb.resize(3);

    Dkb[0] = 1.0 / (e0len * (1.0 + t0.dot(t1)))
      * (2.0 * crossMat(t1) + outerProd(kb, t0 + t1));

    Dkb[2] = 1.0 / (e1len * (1.0 + t0.dot(t1)))
      * (2.0 * crossMat(t0) - outerProd(kb, t0 + t1));

    Dkb[1] = -(Dkb[0] + Dkb[2]);
  }
}

Scalar RodBendingForce::computeDenominator(const Vec3d& x0,
                                           const Vec3d& x1,
                                           const Vec3d& x2)
{
  Vec3d e0 = x1 - x0;
  Vec3d e1 = x2 - x1;
  return e0.norm() * e1.norm() + e0.dot(e1);
}

RodBendingForce::Vec9d RodBendingForce::computeDxDenominator(const Vec3d& x0,
                                                             const Vec3d& x1,
                                                             const Vec3d& x2)
{
  Vec3d e0 = x1 - x0;
  Vec3d e1 = x2 - x1;

  Vec9d DxDenominator;

  DxDenominator.segment(0,3) = -e1 - e1.norm() / e0.norm() * e0;
  DxDenominator.segment(6,3) =  e0 + e0.norm() / e1.norm() * e1;
  DxDenominator.segment(3,3) = -(DxDenominator.segment(0,3)+DxDenominator.segment(6,3));

  return DxDenominator;
}

Scalar RodBendingForce::computeKbDotRef(const Vec3d& x0,
                                        const Vec3d& x1,
                                        const Vec3d& x2,
                                        const Vec3d& u)
{
  Vec3d kb = computeKb(x0, x1, x2);
  return kb.dot(u);
}

RodBendingForce::Vec9d RodBendingForce::computeDxKbDotRef(const Vec3d& x0,
                                                          const Vec3d& x1,
                                                          const Vec3d& x2,
                                                          const Vec3d& u)
{
  Scalar d = computeDenominator(x0, x1, x2);
  Vec9d DxDenominator = computeDxDenominator(x0, x1, x2);
  Vec3d e0 = x1 - x0;
  Vec3d e1 = x2 - x1;

  Vec9d DxKbDotRef;

  DxKbDotRef.segment(0,3) = 2.0*u.cross(e1)/d;
  DxKbDotRef.segment(6,3) = 2.0*u.cross(e0)/d;
  DxKbDotRef.segment(3,3) = -DxKbDotRef.segment(0,3)-DxKbDotRef.segment(6,3);

  DxKbDotRef -= 2.0 * u.dot(e0.cross(e1)) / square(d) * DxDenominator;

  return DxKbDotRef;
}

void RodBendingForce::computeDxDxKbDotRef()
{
  m_DxDxKbDotRef.resize(m_rod.nv());

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    int idx = vh.idx();
    edge_handle eh0 = m_stencil.inEdge();
    edge_handle eh1 = m_stencil.outEdge();
    const Vec3d& e0 = m_rod.getEdge(eh0);
    const Vec3d& e1 = m_rod.getEdge(eh1);
    assert(approxEq(e0, m_rod.getEdge(idx-1)));
    assert(approxEq(e1, m_rod.getEdge(idx)));

    m_DxDxKbDotRef[idx].resize(4);
    const Vec3d& x0 = m_rod.getVertex(idx - 1);
    const Vec3d& x1 = m_rod.getVertex(idx);
    const Vec3d& x2 = m_rod.getVertex(idx + 1);
    Scalar d = m_rod.property(m_denominators)[vh];
    Scalar dSquared = d*d;
    Vec9d DxDenom = computeDxDenominator(x0, x1, x2);
    Mat9d& DxDxDenom = m_DxDxDenominator[idx];

    Mat9d temp = DxDxDenom - 2.0 / d * DxDenom * DxDenom.transpose();
    for (int j = idx - 1; j <= idx; ++j) {
      int j2 = j + 1 - idx;
      const Vec3d& ek = (j2 == 0 ? e0 : e1);
      Scalar div = ek.squaredNorm();

      for (int k = 0; k <= 1; ++k) {
        Mat9d& DxDxKbDotRef = m_DxDxKbDotRef[idx][2 * j2 + k];
        const Vec3d& u = ( k == 0 ?
                           m_rod.getReferenceDirector1(j) :
                           m_rod.getReferenceDirector2(j) );

        DxDxKbDotRef = -2.0 * u.dot(e0.cross(e1)) / dSquared * temp;

        Vec9d myTemp;
        myTemp.segment<3>(0) = -2.0 * u.cross(e1) / dSquared;
        myTemp.segment<3>(6) = -2.0 * u.cross(e0) / dSquared;
        myTemp.segment<3>(3) = -(myTemp.segment<3>(0) + myTemp.segment<3>(6));

        DxDxKbDotRef += DxDenom * myTemp.transpose()
          + myTemp * DxDenom.transpose();

        Mat3d term1 = 2.0 / d * crossMat(u);
        Mat3d term2 = 2.0 / (d * div) * e0.cross(ek) * u.transpose();
        Mat3d term3 = 2.0 / (d * div) * ek.cross(e1) * u.transpose();

        DxDxKbDotRef.block(0,3,3,3) -= term1;
        DxDxKbDotRef.block(0,6,3,3) += term1;
        DxDxKbDotRef.block(3,0,3,3) += term1;
        DxDxKbDotRef.block(3,6,3,3) -= term1;
        DxDxKbDotRef.block(6,0,3,3) -= term1;
        DxDxKbDotRef.block(6,3,3,3) += term1;
        DxDxKbDotRef.block(3,3*j2,3,3) -= term3 - term2;
        DxDxKbDotRef.block(6,3*j2,3,3) -= term2;
        DxDxKbDotRef.block(0,3*j2,3,3) += term3;
        DxDxKbDotRef.block(3,3*(j2+1),3,3) += term3 - term2;
        DxDxKbDotRef.block(6,3*(j2+1),3,3) += term2;
        DxDxKbDotRef.block(0,3*(j2+1),3,3) -= term3;
      }
    }
  }
}

RodBendingForce::Mat9d RodBendingForce::computeDxDxKbDotRef(const Vec3d& x0,
                                                            const Vec3d& x1,
                                                            const Vec3d& x2,
                                                            const Vec3d& u,
                                                            int k,
                                                            int vertIdx)
{
  Scalar d = computeDenominator(x0, x1, x2);
  Scalar dSquared = d*d;
  Vec9d DxDenominator = computeDxDenominator(x0, x1, x2);
  Mat9d& DxDxDenominator = m_DxDxDenominator[vertIdx];
  const Vec3d& e0 = m_rod.getEdge(vertIdx-1);
  const Vec3d& e1 = m_rod.getEdge(vertIdx);

  Mat9d DxDxKbDotRef = -2.0 * u.dot(e0.cross(e1)) / dSquared
    * (DxDxDenominator - 2.0 / d * DxDenominator * DxDenominator.transpose());

  Vec9d myTemp;
  myTemp.segment<3>(0) = -2.0 * u.cross(e1) / dSquared;
  myTemp.segment<3>(6) = -2.0 * u.cross(e0) / dSquared;
  myTemp.segment<3>(3) = -(myTemp.segment<3>(0) + myTemp.segment<3>(6));

  DxDxKbDotRef += DxDenominator * myTemp.transpose() + myTemp * DxDenominator.transpose();

  const Vec3d& ek = (k == 0 ? e0 : e1);
  Scalar div = ek.squaredNorm();
  Mat3d term1 = 2.0 / d * crossMat(u);
  Mat3d term2 = 2.0 / (d * div) * e0.cross(ek) * u.transpose();
  Mat3d term3 = 2.0 / (d * div) * ek.cross(e1) * u.transpose();

  DxDxKbDotRef.block(0,3,3,3) -= term1;
  DxDxKbDotRef.block(0,6,3,3) += term1;
  DxDxKbDotRef.block(3,0,3,3) += term1;
  DxDxKbDotRef.block(3,6,3,3) -= term1;
  DxDxKbDotRef.block(6,0,3,3) -= term1;
  DxDxKbDotRef.block(6,3,3,3) += term1;
  DxDxKbDotRef.block(3,3*k,3,3) -= term3 - term2;
  DxDxKbDotRef.block(6,3*k,3,3) -= term2;
  DxDxKbDotRef.block(0,3*k,3,3) += term3;
  DxDxKbDotRef.block(3,3*(k+1),3,3) += term3 - term2;
  DxDxKbDotRef.block(6,3*(k+1),3,3) += term2;
  DxDxKbDotRef.block(0,3*(k+1),3,3) -= term3;

  return DxDxKbDotRef;
}

Scalar RodBendingForce::globalEnergy()
{
  Scalar energy = 0;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    Scalar localE = localEnergy(vh);
    energy += localE;

#ifdef TEST_ROD_BENDING
    testEnergy(localE, vh);
#endif // TEST_ROD_BENDING
  }

  return energy;
}

Scalar RodBendingForce::localEnergy(const vertex_handle& vh)
{
  edge_handle eh0 = m_stencil.inEdge();
  edge_handle eh1 = m_stencil.outEdge();
  Mat2d B = (getB(eh0) + getB(eh1)) / 2.0;
  Scalar len = getRefVertexLength(vh);

  Scalar e = 0.0;
  for (int j = 0; j <= 1; ++j) {
    const Vec2d& omega = getOmega(vh, j);
    const Vec2d& omegaBar = getOmegaBar(vh, j);
    e += (omega - omegaBar).dot(B * (omega - omegaBar));
  }
  e /= 4.0 * len;

  return e;
}

void RodBendingForce::globalForce(VecXd& force)
{
#ifdef TEST_ROD_BENDING
  globalEnergy();
#endif // TEST_ROD_BENDING

  ElementForce f;
  IndexArray indices;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    f.setZero();
    localForce(f, vh);
    m_stencil.indices(indices);
    for (int j = 0; j < f.size(); ++j) {
      force(indices[j]) += f(j);
    }

#ifdef TEST_ROD_BENDING
    testForce(f, vh);
#endif // TEST_ROD_BENDING
  }
}

void RodBendingForce::localForce(ElementForce& force, const vertex_handle& vh)
{
  edge_handle eh0 = m_stencil.inEdge();
  edge_handle eh1 = m_stencil.outEdge();
  vertex_handle vh0 = m_stencil.prevVert();
  vertex_handle vh1 = m_stencil.nextVert();

  Mat2d B = (getB(eh0) + getB(eh1)) / 2.0;
  //const Vec3d& x0 = m_rod.getVertex(vh0);
  //const Vec3d& x1 = m_rod.getVertex(vh);
  //const Vec3d& x2 = m_rod.getVertex(vh1);
  Scalar len = getRefVertexLength(vh);

  for (int j = 0; j <= 1; ++j) {
    force(4 * j + 3) = -0.5 / len
      * (getOmega(vh, j) - getOmegaBar(vh, j)).dot(B * Jt * getOmega(vh, j));
  }

  int i = vh.idx();
  Eigen::Matrix<Scalar, 1, 3> temp(Eigen::Matrix<Scalar, 1, 3>::Zero().eval());
  for (int j = i - 1; j <= i; ++j) {
    const MaterialMatrix& m = getMaterialMatrix(j);
    int j2 = j + 1 - i;
    temp
      += -0.5 / len * (getOmega(vh, j2) - getOmegaBar(vh, j2)).transpose()
      * B * m;
  }

  const Mat3dArray& Dkb = m_rod.property(m_derivKb)[vh];
  for (int k = 0; k < 3; ++k) {
    force.segment(4 * k, 3) += (temp * Dkb[k]).transpose();
  }
}

void RodBendingForce::globalJacobian(int baseidx, Scalar scale, MatrixBase& Jacobian)
{
  computeDxDenominator();
  computeDxDxDenominator();
  computeDxDxKbDotRef();

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
    Jacobian.add(indices, indices, adder);

#ifdef TEST_ROD_BENDING
    testJacobian(localJ, vh);
#endif // TEST_ROD_BENDING
  }

  /*
#ifdef TIMING_ON
  double comp = Timer::getTimer("computation").getTotal();
  double stuff = Timer::getTimer("stuff").getTotal();
  double loop = Timer::getTimer("loop").getTotal();
  std::cout << stuff/comp << "\t" << loop/comp << std::endl;
#endif // TIMING_ON
  */
}

void RodBendingForce::localJacobian(ElementJacobian& jac,
                                    const vertex_handle& vh)
{
  int idx = vh.idx();
  Mat2d B = (getB(idx - 1) + getB(idx)) / 2.0;
  Scalar len = getRefVertexLength(idx);

  const Mat3dArray& Dkb = m_rod.property(m_derivKb)[vh];

  // \frac{ \partial^2 E }{ \partial^2\theta^j }
  Mat2d JBJt =  J * B * Jt;
  for (int j = idx - 1; j <= idx; ++j) {
    int j2 = j + 1 - idx;
    const Vec2d& omega = getOmega(idx, j);
    const Vec2d& omegaBar = getOmegaBar(idx, j);
    jac(3 + 4 * j2, 3 + 4 * j2) =
      -0.5 / len * (omega.dot(JBJt * omega)
                    - (omega - omegaBar).dot(B * omega));
  }
  jac(3,7) = jac(7,3) = 0;

  Mat9d DxDxKbDotU, DxDxKbDotV;
  for (int j = idx - 1; j <= idx; ++j) {
    int j2 = j + 1 - idx;
    const MaterialMatrix& matVec = getMaterialMatrix(j);
    const Mat2d& Rot = m_rod.property(m_edgeRotationMatrix)[j];
    const Vec2d& omega = getOmega(idx, j);
    const Vec2d& omegaBar = getOmegaBar(idx, j);

    Eigen::Matrix<Scalar, 1, 3> term1
      = -0.5 / len * ((omega - omegaBar).transpose() * B * Jt
                      + omega.transpose() * J * B) * matVec;

    Vec2d term2 = -0.5 / len * (omega - omegaBar).transpose() * B * Rot;
    /*if (term2.isZero()) {
      DxDxKbDotU.setZero();
      DxDxKbDotV.setZero();
      } else {*/
      DxDxKbDotU = m_DxDxKbDotRef[idx][2 * (j + 1 - idx) + 0];
      DxDxKbDotV = m_DxDxKbDotRef[idx][2 * (j + 1 - idx) + 1];
      //}

    for (int i = 0; i < 3; ++i) {
      // \frac{ \partial^2 E }{ \partial\x_i \partial\theta^j }
      jac.block(4 * j2 + 3, 4 * i, 1, 3) = term1 * Dkb[i];
      jac.block(4*i,4*j2+3,3,1) = jac.block(4*j2+3,4*i,1,3).transpose();

      for (int k = 0; k < 3; ++k) {
        // \frac{ \partial^2 E }{ \partial\x_i \partial\x_k }
        jac.block(4*i,4*k,3,3) += term2(0) * DxDxKbDotU.block(3*i,3*k,3,3)
          + term2(1) * DxDxKbDotV.block(3*i,3*k,3,3);
      }
    }
  }

  Mat3d accum(Mat3d::Zero());
  for (int j = idx - 1; j <= idx; ++j) {
    const MaterialMatrix& matVec = getMaterialMatrix(j);
    accum += matVec.transpose() * B * matVec;
  }
  accum *= -0.5 / len;
  for (int i = 0; i < 3; ++i) {
    for (int k = 0; k < 3; ++k) {
      // \frac{ \partial^2 E }{ \partial\x_i \partial\x_k }
      jac.block(4*i,4*k,3,3) += Dkb[i].transpose() * accum * Dkb[k];
    }
  }
}

#ifdef TEST_ROD_BENDING

void RodBendingForce::testEnergy(const Scalar& energy,
                                 const vertex_handle& vh) const
{
  const Vec3d& x0 = m_rod.getVertex(m_stencil.prevVert());
  const Vec3d& x1 = m_rod.getVertex(vh);
  const Vec3d& x2 = m_rod.getVertex(m_stencil.nextVert());
  const Vec2dArray& omegaBar = m_rod.property(m_omegaBar)[vh];
  Scalar theta0 = m_rod.getTheta(m_stencil.inEdge());
  Scalar theta1 = m_rod.getTheta(m_stencil.outEdge());
  const Vec3d& u0 = m_rod.getReferenceDirector1(m_stencil.inEdge());
  const Vec3d& u1 = m_rod.getReferenceDirector1(m_stencil.outEdge());
  const Vec3d& v0 = m_rod.getReferenceDirector2(m_stencil.inEdge());
  const Vec3d& v1 = m_rod.getReferenceDirector2(m_stencil.outEdge());
  Mat2d B = (getB(m_stencil.inEdge()) + getB(m_stencil.outEdge())) / 2.0;
  Scalar len = getRefVertexLength(vh);
  Scalar mathE;
  rodBendingEnergyTest(mathE, energy, x0, x1, x2, theta0, theta1,
                       omegaBar[0], omegaBar[1], u0, u1, v0, v1, B, len);
}

void RodBendingForce::testForce(const ElementForce& force,
                                const vertex_handle& vh) const
{
  const Vec3d& x0 = m_rod.getVertex(m_stencil.prevVert());
  const Vec3d& x1 = m_rod.getVertex(vh);
  const Vec3d& x2 = m_rod.getVertex(m_stencil.nextVert());
  Scalar theta0 = m_rod.getTheta(m_stencil.inEdge());
  Scalar theta1 = m_rod.getTheta(m_stencil.outEdge());
  const Vec3d& u0 = m_rod.getReferenceDirector1(m_stencil.inEdge());
  const Vec3d& u1 = m_rod.getReferenceDirector1(m_stencil.outEdge());
  const Vec3d& v0 = m_rod.getReferenceDirector2(m_stencil.inEdge());
  const Vec3d& v1 = m_rod.getReferenceDirector2(m_stencil.outEdge());
  Mat2d B = (getB(m_stencil.inEdge()) + getB(m_stencil.outEdge())) / 2.0;
  const Vec2dArray& omegaBar = m_rod.property(m_omegaBar)[vh];
  Scalar len = getRefVertexLength(vh);
  ElementForce mathF;
  if (rodBendingForceTest(mathF, force, x0, x1, x2, theta0, theta1,
                          omegaBar[0], omegaBar[1], u0, u1, v0, v1, B, len)) {
    std::cout << "mathF = " << mathF << std::endl;
    std::cout << "myF = " << force << std::endl << std::endl;
  }
}

void RodBendingForce::testJacobian(const ElementJacobian& Jacobian,
                                   const vertex_handle& vh) const
{
  const Vec3d& x0 = m_rod.getVertex(m_stencil.prevVert());
  const Vec3d& x1 = m_rod.getVertex(vh);
  const Vec3d& x2 = m_rod.getVertex(m_stencil.nextVert());
  Scalar theta0 = m_rod.getTheta(m_stencil.inEdge());
  Scalar theta1 = m_rod.getTheta(m_stencil.outEdge());
  const Vec3d& u0 = m_rod.getReferenceDirector1(m_stencil.inEdge());
  const Vec3d& u1 = m_rod.getReferenceDirector1(m_stencil.outEdge());
  const Vec3d& v0 = m_rod.getReferenceDirector2(m_stencil.inEdge());
  const Vec3d& v1 = m_rod.getReferenceDirector2(m_stencil.outEdge());
  Mat2d B = (getB(m_stencil.inEdge()) + getB(m_stencil.outEdge())) / 2.0;
  const Vec2dArray& omegaBar = m_rod.property(m_omegaBar)[vh];
  Scalar len = getRefVertexLength(vh);
  ElementJacobian mathJ;
  if (rodBendingJacobianTest(mathJ, Jacobian, x0, x1, x2, theta0, theta1,
                             omegaBar[0], omegaBar[1], u0, u1, v0, v1,
                             B, len) == -1) {
    std::cout << "localJ = " << Jacobian << ";" << std::endl;
    std::cout << "mathJ = " << mathJ << ";" << std::endl;
  }
}

#endif // TEST_ROD_BENDING

} // namespace BASim

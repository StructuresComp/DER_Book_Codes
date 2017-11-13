/**
 * \file RodForce.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/01/2009
 */

#include "BASim/src/Physics/ElasticRods/RodForce.hh"
#include "BASim/src/Math/Math.hh"

namespace BASim {

using namespace Util;

Mat2d RodForce::J(Mat2d::Zero());
Mat2d RodForce::Jt(Mat2d::Zero());

RodForce::RodForce(ElasticRod& rod, const std::string& name)
  : m_rod(rod)
  , m_name(name)
  , m_viscous(false)
{
  J << 0, -1, 1, 0;
  Jt << 0, 1, -1, 0;
}

const std::string& RodForce::getName() const
{
  return m_name;
}

void RodForce::computeKb(Vec3d& kb, const Vec3d& x0, const Vec3d& x1,
                         const Vec3d& x2)
{
  Vec3d e0 = x1 - x0;
  Vec3d e1 = x2 - x1;
  kb = 2.0 * e0.cross(e1) / (e0.norm() * e1.norm() + e0.dot(e1));
}

Vec3d RodForce::computeKb(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2)
{
  Vec3d kb;
  computeKb(kb, x0, x1, x2);
  return kb;
}

void RodForce::computeDkb(Mat3dArray& Dkb, const Vec3d& x0, const Vec3d& x1,
                          const Vec3d& x2)
{
  Dkb.resize(3);

  Vec3d e0 = x1 - x0;
  Vec3d t0 = e0.normalized();
  Vec3d e1 = x2 - x1;
  Vec3d t1 = e1.normalized();

  Vec3d kb = computeKb(x0, x1, x2);

  Dkb[0] = 1.0/(e0.norm() * (1.0 + t0.dot(t1)))
    * (2.0 * crossMat(t1) + outerProd(kb, t0 + t1));

  Dkb[2] = 1.0/(e1.norm() * (1.0 + t0.dot(t1)))
    * (2.0 * crossMat(t0) - outerProd(kb, t0 + t1));

  Dkb[1] = -(Dkb[0] + Dkb[2]);

#ifdef VERIFY_FORCES
  Scalar xi = x0(0);
  Scalar yi = x0(1);
  Scalar zi = x0(2);
  Scalar xj = x1(0);
  Scalar yj = x1(1);
  Scalar zj = x1(2);
  Scalar xk = x2(0);
  Scalar yk = x2(1);
  Scalar zk = x2(2);
  MatXd test(3,9);
#include "derivKb.inl"
  for (int i = 0; i < 3; ++i) {
    if((Dkb[i]-test.block(0,3*i,3,3)).norm() > SMALL_NUMBER) {
      std::cout << "ERROR IN DKB" << std::endl;
      std::cout << Dkb[i] << std::endl
                << test.block(0,3*i,3,3) << std::endl;
    }
  }
#endif // VERIFY_FORCES
}

Mat3dArray RodForce::computeDkb(const Vec3d& x0, const Vec3d& x1,
                                const Vec3d& x2)
{
  Mat3dArray Dkb;
  computeDkb(Dkb, x0, x1, x2);
  return Dkb;
}

} // namespace BASim

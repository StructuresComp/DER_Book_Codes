/**
 * \file RodBendingTest.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/07/2009
 */

#include "RodBendingTest.hh"
#include "BASim/src/Core/Util.hh"
#include "mdefs.h"

namespace BASim {

int rodBendingEnergyTest(Scalar& mathematicaEnergy,
                         const Scalar& energy,
                         const Vec3d& x0,
                         const Vec3d& x1,
                         const Vec3d& x2,
                         const Scalar& theta0,
                         const Scalar& theta1,
                         const Vec2d& omegaBar0,
                         const Vec2d& omegaBar1,
                         const Vec3d& u0,
                         const Vec3d& u1,
                         const Vec3d& v0,
                         const Vec3d& v1,
                         const Mat2d& B,
                         const Scalar& len)
{
  Scalar xi = x0(0);
  Scalar yi = x0(1);
  Scalar zi = x0(2);
  Scalar xj = x1(0);
  Scalar yj = x1(1);
  Scalar zj = x1(2);
  Scalar xk = x2(0);
  Scalar yk = x2(1);
  Scalar zk = x2(2);
  Scalar omegaBar11 = omegaBar0(0);
  Scalar omegaBar12 = omegaBar0(1);
  Scalar omegaBar21 = omegaBar1(0);
  Scalar omegaBar22 = omegaBar1(1);
  Scalar ux0 = u0(0);
  Scalar uy0 = u0(1);
  Scalar uz0 = u0(2);
  Scalar ux1 = u1(0);
  Scalar uy1 = u1(1);
  Scalar uz1 = u1(2);
  Scalar vx0 = v0(0);
  Scalar vy0 = v0(1);
  Scalar vz0 = v0(2);
  Scalar vx1 = v1(0);
  Scalar vy1 = v1(1);
  Scalar vz1 = v1(2);
  Scalar kb1 = B(0,0);
  Scalar kb2 = B(1,1);
#include "bendingEnergy.inl"
  if (!approxEq(energy, mathematicaEnergy)) {
    std::cerr << "Error in bending energy" << std::endl
              << energy << "\t" << mathematicaEnergy << std::endl;
    return -1;
  }
  return 0;
}

int rodBendingForceTest(Eigen::Matrix<Scalar, 11, 1>& mathematicaForce,
                        const Eigen::Matrix<Scalar, 11, 1>& force,
                        const Vec3d& x0,
                        const Vec3d& x1,
                        const Vec3d& x2,
                        const Scalar& theta0,
                        const Scalar& theta1,
                        const Vec2d& omegaBar0,
                        const Vec2d& omegaBar1,
                        const Vec3d& u0,
                        const Vec3d& u1,
                        const Vec3d& v0,
                        const Vec3d& v1,
                        const Mat2d& B,
                        const Scalar& len)
{
  Scalar xi = x0(0);
  Scalar yi = x0(1);
  Scalar zi = x0(2);
  Scalar xj = x1(0);
  Scalar yj = x1(1);
  Scalar zj = x1(2);
  Scalar xk = x2(0);
  Scalar yk = x2(1);
  Scalar zk = x2(2);
  Scalar omegaBar11 = omegaBar0(0);
  Scalar omegaBar12 = omegaBar0(1);
  Scalar omegaBar21 = omegaBar1(0);
  Scalar omegaBar22 = omegaBar1(1);
  Scalar ux0 = u0(0);
  Scalar uy0 = u0(1);
  Scalar uz0 = u0(2);
  Scalar ux1 = u1(0);
  Scalar uy1 = u1(1);
  Scalar uz1 = u1(2);
  Scalar vx0 = v0(0);
  Scalar vy0 = v0(1);
  Scalar vz0 = v0(2);
  Scalar vx1 = v1(0);
  Scalar vy1 = v1(1);
  Scalar vz1 = v1(2);
  Scalar kb1 = B(0,0);
  Scalar kb2 = B(1,1);
#include "bendingForce.inl"
  if (!approxEq(force, mathematicaForce)) {
    std::cerr << "Error in bending force" << std::endl;
    return -1;
  }
  return 0;
}

int rodBendingJacobianTest(Eigen::Matrix<Scalar, 11, 11>& mathematicaJacobian,
                           const Eigen::Matrix<Scalar, 11, 11>& Jacobian,
                           const Vec3d& x0,
                           const Vec3d& x1,
                           const Vec3d& x2,
                           const Scalar& theta0,
                           const Scalar& theta1,
                           const Vec2d& omegaBar0,
                           const Vec2d& omegaBar1,
                           const Vec3d& u0,
                           const Vec3d& u1,
                           const Vec3d& v0,
                           const Vec3d& v1,
                           const Mat2d& B,
                           const Scalar& len)
{
  Scalar xi = x0(0);
  Scalar yi = x0(1);
  Scalar zi = x0(2);
  Scalar xj = x1(0);
  Scalar yj = x1(1);
  Scalar zj = x1(2);
  Scalar xk = x2(0);
  Scalar yk = x2(1);
  Scalar zk = x2(2);
  Scalar omegaBar11 = omegaBar0(0);
  Scalar omegaBar12 = omegaBar0(1);
  Scalar omegaBar21 = omegaBar1(0);
  Scalar omegaBar22 = omegaBar1(1);
  Scalar ux0 = u0(0);
  Scalar uy0 = u0(1);
  Scalar uz0 = u0(2);
  Scalar ux1 = u1(0);
  Scalar uy1 = u1(1);
  Scalar uz1 = u1(2);
  Scalar vx0 = v0(0);
  Scalar vy0 = v0(1);
  Scalar vz0 = v0(2);
  Scalar vx1 = v1(0);
  Scalar vy1 = v1(1);
  Scalar vz1 = v1(2);
  Scalar kb1 = B(0,0);
  Scalar kb2 = B(1,1);
#include "bendingJacobian.inl"
  if (!approxEq(Jacobian, mathematicaJacobian)) {
    std::cerr << "Error in bending Jacobian" << std::endl;
    return -1;
  }
  return 0;
}

} // namespace BASim

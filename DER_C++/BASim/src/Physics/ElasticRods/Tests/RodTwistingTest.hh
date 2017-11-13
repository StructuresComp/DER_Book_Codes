/**
 * \file RodTwistingTest.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/07/2009
 */

#ifndef RODTWISTINGTEST_HH
#define RODTWISTINGTEST_HH

#include "mdefs.h"

namespace BASim {

inline int rodTwistingEnergyTest(Scalar& mathematicaEnergy,
                                 const Scalar& energy,
                                 const Scalar& phi,
                                 const Scalar& theta0,
                                 const Scalar& theta1,
                                 const Scalar& undefTwist,
                                 const Scalar& kt,
                                 const Scalar& len)
{
#include "twistingEnergy.inl"
  if (!approxEq(energy, mathematicaEnergy)) {
    std::cerr << "Error in twisting energy" << std::endl;
    return -1;
  }
  return 0;
}

inline int rodTwistingForceTest(Eigen::Matrix<Scalar, 11, 1>& mathematicaForce,
                                const Eigen::Matrix<Scalar, 11, 1>& force,
                                const Vec3d& x0,
                                const Vec3d& x1,
                                const Vec3d& x2,
                                const Scalar& theta0,
                                const Scalar& theta1,
                                const Scalar& phi,
                                const Scalar& undefTwist,
                                const Scalar& kt,
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
#include "twistingForce.inl"
  if (!approxEq(force, mathematicaForce)) {
    std::cerr << "Error in twisting force" << std::endl;
    return -1;
  }
  return 0;
}

inline int rodTwistingJacobianTest(Eigen::Matrix<Scalar, 11, 11>& mathematicaJacobian,
                                   const Eigen::Matrix<Scalar, 11, 11>& Jacobian,
                                   const Vec3d& x0,
                                   const Vec3d& x1,
                                   const Vec3d& x2,
                                   const Scalar& theta0,
                                   const Scalar& theta1,
                                   const Scalar& phi,
                                   const Scalar& undefTwist,
                                   const Scalar& kt,
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
#include "twistingJacobian.inl"
  if (!approxEq(Jacobian, mathematicaJacobian)) {
    std::cerr << "Error in twisting Jacobian" << std::endl;
    return -1;
  }
  return 0;
}

} // namespace BASim

#endif // RODTWISTINGTEST_HH

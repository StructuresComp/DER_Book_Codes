/**
 * \file RodStretchingTest.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/07/2009
 */

#ifndef RODSTRETCHINGTEST_HH
#define RODSTRETCHINGTEST_HH

#include "mdefs.h"

namespace BASim {

inline int
rodStretchingEnergyTest(Scalar& mathematicaEnergy,
                        const Scalar& energy,
                        const Vec3d& x0,
                        const Vec3d& x1,
                        const Scalar& ks,
                        const Scalar& len)
{
  Scalar xi = x0(0);
  Scalar yi = x0(1);
  Scalar zi = x0(2);
  Scalar xj = x1(0);
  Scalar yj = x1(1);
  Scalar zj = x1(2);
#include "stretchingEnergy.inl"
  if (!approxEq(energy, mathematicaEnergy)) {
    std::cerr << "Error in stretching energy" << std::endl;
    return -1;
  }
  return 0;
}

inline int
rodStretchingForceTest(Eigen::Matrix<Scalar, 6, 1>& mathematicaForce,
                       const Eigen::Matrix<Scalar, 6, 1>& force,
                       const Vec3d& x0,
                       const Vec3d& x1,
                       const Scalar& ks,
                       const Scalar& len)
{
  Scalar xi = x0(0);
  Scalar yi = x0(1);
  Scalar zi = x0(2);
  Scalar xj = x1(0);
  Scalar yj = x1(1);
  Scalar zj = x1(2);
#include "stretchingForce.inl"
  if (!approxEq(force, mathematicaForce)) {
    std::cerr << "Error in stretching force" << std::endl;
    return -1;
  }
  return 0;
}

inline int
rodStretchingJacobianTest(Eigen::Matrix<Scalar, 6, 6>& mathematicaJacobian,
                          const Eigen::Matrix<Scalar, 6, 6>& Jacobian,
                          const Vec3d& x0,
                          const Vec3d& x1,
                          const Scalar& ks,
                          const Scalar& len)
{
  Scalar xi = x0(0);
  Scalar yi = x0(1);
  Scalar zi = x0(2);
  Scalar xj = x1(0);
  Scalar yj = x1(1);
  Scalar zj = x1(2);
#include "stretchingJacobian.inl"
  if (!approxEq(Jacobian, mathematicaJacobian)) {
    std::cerr << "Error in stretching Jacobian" << std::endl;
    return -1;
  }
  return 0;
}

} // namespace BASim

#endif // RODSTRETCHINGTEST_HH

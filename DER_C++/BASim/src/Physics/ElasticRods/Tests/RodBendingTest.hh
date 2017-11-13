/**
 * \file RodBendingTest.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/07/2009
 */

#ifndef RODBENDINGTEST_HH
#define RODBENDINGTEST_HH

#include "BASim/src/Core/Definitions.hh"

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
                         const Scalar& len);

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
                        const Scalar& len);

int rodBendingJacobianTest(Eigen::Matrix<Scalar, 11, 11>& mathematicaJacobian,
                           const Eigen::Matrix<Scalar, 11, 11>& jacobian,
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
                           const Scalar& len);

} // namespace BASim

#endif // RODBENDINGTEST_HH

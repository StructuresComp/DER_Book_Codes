/**
 * \file RodAnisoForce.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/11/2009
 */

#ifndef RODANISOFORCE_HH
#define RODANISOFORCE_HH

#include "BASim/src/Physics/ElasticRods/RodAnisoBending.hh"

namespace BASim {

/** Class that computes forces for the rod assuming the reference
    frame is given by the Bishop (space parallel) frame. */
class RodAnisoForce : public RodAnisoBending
{
public:
  RodAnisoForce(ElasticRod& rod);

  Scalar globalEnergy();
  void globalForce(VecXd& force);
  void globalJacobian(int baseidx, Scalar scale, MatrixBase& Jacobian);
  /*
  void grad(VecXd& G);
  void hess(Matrix& H);

  void checkGrad();
  void checkHess();

  void evalTorque(Vec3d& t, int j);
  */

  void updateProperties();

protected:
  void computeNablaKappa();
  void evalNablaKappa(Mat3d& nk, int i, int v);
  Mat3d& evalNablaKappa(int i, int v);
  void fdNablaKappa(Mat3d& nk, int i, int v);
  void computeNablaPsi();
  void evalNablaPsi(Vec3d& np, int i, int v);
  Vec3d& evalNablaPsi(int i, int v);
  void fdNablaPsi(Vec3d& np, int i, int v);

  std::vector<Mat3d> _nk;
  std::vector<Vec3d> _np;

  VecXd forPsi;
};

} // namespace BASim

#endif // RODANISOFORCE_HH

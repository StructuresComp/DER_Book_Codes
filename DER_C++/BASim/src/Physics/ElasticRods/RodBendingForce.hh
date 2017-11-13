/**
 * \file RodBendingForce.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/06/2009
 */

#ifndef RODBENDINGFORCE_HH
#define RODBENDINGFORCE_HH

#include "BASim/src/Physics/ElasticRods/RodAnisoBending.hh"

namespace BASim {

/** This class implements the anisotropic bending force for an elastic
    rod. */
class RodBendingForce : public RodAnisoBending
{
public:

  typedef Eigen::Matrix<Scalar, 9, 1> Vec9d;
  typedef Eigen::Matrix<Scalar, 9, 9> Mat9d;

  RodBendingForce(ElasticRod& rod);

  virtual Scalar globalEnergy();
  virtual void globalForce(VecXd& force);
  virtual void globalJacobian(int baseidx, Scalar scale, MatrixBase& J);

  Scalar localEnergy(const vertex_handle& vh);
  void localForce(ElementForce& force, const vertex_handle& vh);
  void localJacobian(ElementJacobian& J, const vertex_handle& vh);

  Scalar computeDenominator(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2);
  Vec9d computeDxDenominator(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2);
  //Mat9d computeDxDxDenominator(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2);

  Scalar computeKbDotRef(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2, const Vec3d& u);
  Vec9d computeDxKbDotRef(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2, const Vec3d& u);
  Mat9d computeDxDxKbDotRef(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2, const Vec3d& u, int k, int vertIdx);
  void computeDxDxKbDotRef();

  virtual void updateProperties();
  void computeRotationMatrix();
  void computeDenominators();
  void computeDxDenominator();
  void computeDxDxDenominator();

  void computeDkb();

  using RodAnisoBending::computeDkb;

protected:

#ifdef TEST_ROD_BENDING
  void testEnergy(const Scalar& energy, const vertex_handle& vh) const;
  void testForce(const ElementForce& force, const vertex_handle& vh) const;
  void testJacobian(const ElementJacobian& Jacobian,
                    const vertex_handle& vh) const;
#endif // TEST_ROD_BENDING

  bool firstDerivs;
  bool secondDerivs;
  std::vector<Vec2dArray> m_kbDotRef;

  VPropHandle<Scalar> m_denominators; ///< denominator of expression for curvature binormal
  std::vector<Vec9d> m_DxDenominator;
  std::vector<Mat9d> m_DxDxDenominator;
  std::vector< std::vector<Mat9d> > m_DxDxKbDotRef;
  VPropHandle<Mat3dArray> m_derivKb; ///< derivative of curvature binormal w.r.t. vertices
  EPropHandle<Mat2d> m_edgeRotationMatrix;
};

} // namespace BASim

#endif // RODBENDINGFORCE_HH

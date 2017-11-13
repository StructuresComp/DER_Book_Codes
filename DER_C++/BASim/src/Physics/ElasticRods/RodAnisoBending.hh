/**
 * \file RodAnisoBending.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/11/2009
 */

#ifndef RODANISOBENDING_HH
#define RODANISOBENDING_HH

#include "BASim/src/Physics/ElasticRods/RodForce.hh"
#include "BASim/src/Physics/ElasticRods/VertexStencil.hh"

namespace BASim {

/** Base class for anisotropic/curved rod bending force. */
class RodAnisoBending : public RodForceT<VertexStencil>
{
public:

  typedef Eigen::Matrix<Scalar, 11, 1> ElementForce;
  typedef Eigen::Matrix<Scalar, 11, 11> ElementJacobian;
  typedef Eigen::Matrix<Scalar, 2, 3> MaterialMatrix;

  RodAnisoBending(ElasticRod& rod, const std::string& name = "RodAnisoBending");

  void computeCurvatureBinormals();
  void computeOmega();
  void computeMaterialMatrix();

  const Mat2d& getB(const edge_handle& eh) const;
  void setB(const edge_handle& eh, const Mat2d& B);

  const Vec3d& getCurvatureBinormal(const vertex_handle& vh) const;
  void setCurvatureBinormal(const vertex_handle& vh, const Vec3d& kb);

  /** returns omega (material curvature vector) given the vertex
      handle and a 0 or 1 indicating the incoming or outgoing edge */
  const Vec2d& getOmega(const vertex_handle& vh, int j) const;
  void setOmega(const vertex_handle& vh, int j, const Vec2d& omega);

  const Vec2d& getOmegaBar(const vertex_handle& vh, int j) const;
  void setOmegaBar(const vertex_handle& vh, int j, const Vec2d& omegaBar);

  Scalar getRefVertexLength(const vertex_handle& vh) const;
  void setRefVertexLength(const vertex_handle& vh, Scalar length);

  const MaterialMatrix& getMaterialMatrix(const edge_handle& eh) const;
  void setMaterialMatrix(const edge_handle& eh, const MaterialMatrix& matrix);
  //--------------------------------------------------------------------------
  const Mat2d& getB(int j) const;
  void setB(int j, const Mat2d& B);

  const Vec3d& getCurvatureBinormal(int i) const;
  void setCurvatureBinormal(int i, const Vec3d& kb);

  const Vec2d& getOmega(int i, int j) const;
  void setOmega(int i, int j, const Vec2d& omega);

  const Vec2d& getOmegaBar(int i, int j) const;
  void setOmegaBar(int i, int j, const Vec2d& omegaBar);

  Scalar getRefVertexLength(int i) const;
  void setRefVertexLength(int i, Scalar length);

  const MaterialMatrix& getMaterialMatrix(int j) const;
  void setMaterialMatrix(int j, const MaterialMatrix& matrix);

  virtual void updateProperties();
  virtual void updateStiffness();
  virtual void updateUndeformedStrain();
  virtual void updateReferenceDomain();

protected:

  EPropHandle<Mat2d> m_B;
  VPropHandle<Scalar> m_refVertexLength;

  VPropHandle<Vec3d> m_curvatureBinormal;
  VPropHandle<Vec2dArray> m_omega;
  VPropHandle<Vec2dArray> m_omegaBar;

  EPropHandle<MaterialMatrix> m_materialMatrix;
};

} // namespace BASim

#endif // RODANISOBENDING_HH

/**
 * \file RodAnisoBending.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/11/2009
 */

#include "BASim/src/Physics/ElasticRods/RodAnisoBending.hh"
#include "BASim/src/Math/Math.hh"

namespace BASim {

RodAnisoBending::RodAnisoBending(ElasticRod& rod, const std::string& name)
  : RodForceT<VertexStencil>(rod, name)
{
  m_rod.add_property(m_B, "bending matrix");
  m_rod.add_property(m_refVertexLength, "bending ref length");

  m_rod.add_property(m_curvatureBinormal, "curvature binormal");
  m_rod.add_property(m_omega, "omega", Vec2dArray(2, Vec2d::Zero()));
  m_rod.add_property(m_omegaBar, "omega bar", Vec2dArray(2, Vec2d::Zero()));
  m_rod.add_property(m_materialMatrix, "material matrix",
                     Eigen::Matrix<Scalar, 2, 3>::Zero().eval());

  updateProperties();
  updateUndeformedStrain();
  updateStiffness();
  updateReferenceDomain();
}

void RodAnisoBending::updateProperties()
{
  computeCurvatureBinormals();
  computeOmega();
  computeMaterialMatrix();
}

void RodAnisoBending::updateUndeformedStrain()
{
  for (int i = 1; i < m_rod.nv() - 1; ++i) {
    for (int j = i - 1; j <= i; ++j) {
      setOmegaBar(i, j, getOmega(i, j));
    }
  }
}

const Mat2d& RodAnisoBending::getB(const edge_handle& eh) const
{
  return m_rod.property(m_B)[eh];
}

void RodAnisoBending::setB(const edge_handle& eh, const Mat2d& B)
{
  m_rod.property(m_B)[eh] = B;
}

Scalar RodAnisoBending::getRefVertexLength(const vertex_handle& vh) const
{
  return m_rod.property(m_refVertexLength)[vh];
}

void RodAnisoBending::setRefVertexLength(const vertex_handle& vh, Scalar length)
{
  m_rod.property(m_refVertexLength)[vh] = length;
}

const Mat2d& RodAnisoBending::getB(int j) const
{
  return m_rod.property(m_B)[j];
}

void RodAnisoBending::setB(int j, const Mat2d& B)
{
  m_rod.property(m_B)[j] = B;
}

const Vec2d& RodAnisoBending::getOmega(const vertex_handle& vh, int j) const
{
  assert(j >= 0);
  assert(j <= 1);
  return m_rod.property(m_omega)[vh][j];
}

void RodAnisoBending::setOmega(const vertex_handle& vh, int j,
                               const Vec2d& omega)
{
  assert(j >= 0);
  assert(j <= 1);
  m_rod.property(m_omega)[vh][j] = omega;
}

const Vec2d& RodAnisoBending::getOmegaBar(const vertex_handle& vh, int j) const
{
  assert(j >= 0);
  assert(j <= 1);
  return m_rod.property(m_omegaBar)[vh][j];
}

void RodAnisoBending::setOmegaBar(const vertex_handle& vh, int j,
                                  const Vec2d& omegaBar)
{
  assert(j >= 0);
  assert(j <= 1);
  m_rod.property(m_omegaBar)[vh][j] = omegaBar;
}

const Vec2d& RodAnisoBending::getOmega(int i, int j) const
{
  return m_rod.property(m_omega)[i][j + 1 - i];
}

void RodAnisoBending::setOmega(int i, int j, const Vec2d& omega)
{
  m_rod.property(m_omega)[i][j + 1 - i] = omega;
}

const Vec2d& RodAnisoBending::getOmegaBar(int i, int j) const
{
  return m_rod.property(m_omegaBar)[i][j + 1 - i];
}

void RodAnisoBending::setOmegaBar(int i, int j, const Vec2d& omegaBar)
{
  m_rod.property(m_omegaBar)[i][j + 1 - i] = omegaBar;
}

Scalar RodAnisoBending::getRefVertexLength(int i) const
{
  return m_rod.property(m_refVertexLength)[i];
}

void RodAnisoBending::setRefVertexLength(int i, Scalar length)
{
  m_rod.property(m_refVertexLength)[i] = length;
}

const Vec3d& RodAnisoBending::getCurvatureBinormal(int i) const
{
  return m_rod.property(m_curvatureBinormal)[i];
}

void RodAnisoBending::setCurvatureBinormal(int i, const Vec3d& kb)
{
  m_rod.property(m_curvatureBinormal)[i] = kb;
}

const Eigen::Matrix<Scalar, 2, 3>&
RodAnisoBending::getMaterialMatrix(int j) const
{
  return m_rod.property(m_materialMatrix)[j];
}

void RodAnisoBending::setMaterialMatrix(int j,
  const Eigen::Matrix<Scalar, 2, 3>& matrix)
{
  m_rod.property(m_materialMatrix)[j] = matrix;
}

void RodAnisoBending::updateStiffness()
{
  for (int j = 0; j < m_rod.ne(); ++j) {
    Scalar a = m_rod.radiusA(j);
    Scalar b = m_rod.radiusB(j);
    Scalar E = m_rod.getYoungsModulus();
    m_rod.property(m_B)[j] << M_PI * pow(a, 4) * E / 4.0, 0, 0,
      M_PI * pow(b, 4) * E / 4.0;
  }
}

void RodAnisoBending::updateReferenceDomain()
{
  for (int i = 1; i < m_rod.nv() - 1; ++i) {
    Scalar len = (m_rod.getEdgeLength(i - 1) + m_rod.getEdgeLength(i)) / 2.0;
    setRefVertexLength(i, len);
  }
}

void RodAnisoBending::computeCurvatureBinormals()
{
  for (int i = 1; i < m_rod.nv() - 1; ++i) {
    Vec3d& kb = m_rod.property(m_curvatureBinormal)[i];
    computeCurvatureBinormal(kb, m_rod.getTangent(i - 1), m_rod.getTangent(i));
  }
}

void RodAnisoBending::computeOmega()
{
  for (int i = 1; i < m_rod.nv() - 1; ++i) {
    const Vec3d& kb = getCurvatureBinormal(i);
    for (int j = i - 1; j <= i; ++j) {
      const Vec3d& m1 = m_rod.getMaterial1(j);
      const Vec3d& m2 = m_rod.getMaterial2(j);
      Vec2d omega(kb.dot(m2), -kb.dot(m1));
      setOmega(i, j, omega);
    }
  }
}

void RodAnisoBending::computeMaterialMatrix()
{
  for (int j = 0; j < m_rod.ne(); ++j) {
    const Vec3d& m1 = m_rod.getMaterial1(j);
    const Vec3d& m2 = m_rod.getMaterial2(j);
    m_rod.property(m_materialMatrix)[j]
      << m2(0),  m2(1),  m2(2), -m1(0), -m1(1), -m1(2);
  }
}

} // namespace BASim

/**
 * \file RodForce.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/01/2009
 */

#ifndef RODFORCE_HH
#define RODFORCE_HH

#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Math/MatrixBase.hh"

namespace BASim {

/** Base class for a force that acts on rods. */
class RodForce
{
public:

  typedef ElasticRod::vertex_handle vertex_handle;
  typedef ElasticRod::edge_handle   edge_handle;

  RodForce(ElasticRod& rod, const std::string& name = "RodForce");
  virtual ~RodForce() {}

  const std::string& getName() const;

  void computeKb(Vec3d& kb, const Vec3d& x0, const Vec3d& x1, const Vec3d& x2);
  Vec3d computeKb(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2);

  void computeDkb(Mat3dArray& Dkb, const Vec3d& x0, const Vec3d& x1,
                  const Vec3d& x2);
  Mat3dArray computeDkb(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2);

  virtual Scalar globalEnergy() = 0;
  virtual void globalForce(VecXd& force) = 0;
  virtual void globalJacobian(int baseidx, Scalar scale, MatrixBase& Jacobian) = 0;

  virtual void updateProperties() {}
  virtual void updateStiffness() {}
  virtual void updateUndeformedStrain() {}
  virtual void updateReferenceDomain() {}

  virtual void verifyProperties() {}
  
  virtual void updatePropertiesForNewVertex() {}

  bool viscous() const { return m_viscous; }
  void setViscous(bool v) { m_viscous = v; updateStiffness(); }

protected:

  ElasticRod& m_rod;
  std::string m_name;
  bool m_viscous;

  static Mat2d J;
  static Mat2d Jt;
};

template <class Stencil>
class RodForceT : public RodForce
{
public:

  typedef typename Stencil::iterator iterator;

  RodForceT(ElasticRod& rod, const std::string& name = "RodForce")
    : RodForce(rod,name)
    , m_stencil(rod)
  {}

protected:

  Stencil m_stencil;
};

} // namespace BASim

#endif // RODFORCE_HH

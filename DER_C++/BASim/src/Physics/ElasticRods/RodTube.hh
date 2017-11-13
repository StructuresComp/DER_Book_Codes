/**
 * \file RodTube.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/11/2009
 */

#ifndef RODTUBE_HH
#define RODTUBE_HH

#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"

namespace BASim {

/** Creates an extruded surface for a rod. */
class RodTube
{
public:

  RodTube(ElasticRod& rod, int slices = 12, bool twistPoints = true);
  ~RodTube() {}

  void setSlices(int slices);
  int getSlices() const;

  bool getTwistPoints() const;
  void setTwistPoints(bool twistPoints);

  void buildTube();
  const std::vector<Vec3d>& getTubeAtVert(int vertIdx) const;
  const Vec3d& getPointAtVert(int vertIdx, int pointIdx) const;
  std::pair<Scalar,Scalar> getUVforPointAtVert(int vertIdx, int pointIdx);

protected:

  void computeAngles();
  void buildFrame();

  ElasticRod& m_rod;
  int m_slices;       ///< number of points to sample cross section with
  bool m_twistPoints; ///< twist the points or only their uv coordinates

  ScalarArray m_sin2; ///< precomputed sine squared
  ScalarArray m_cos2; ///< precomputed cosine square

  EPropHandle<Vec3d> m_frame; ///< material frame at each edge for rendering
  VPropHandle< std::vector<Vec3d> > m_points; ///< extruded points at a vertex
  VPropHandle< std::vector< std::pair<Scalar, Scalar> > > m_uvCoord; ///< texture coordinates at a vertex
  VPropHandle<Scalar> m_twist; ///< twist at a vertex
};

} // namespace BASim

#endif // RODTUBE_HH

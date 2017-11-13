#include "BASim/src/Physics/ElasticRods/RodTube.hh"
#include "BASim/src/Math/Math.hh"

using namespace std;

namespace BASim {

RodTube::RodTube(ElasticRod& rod, int slices, bool twistPoints)
  : m_rod(rod)
  , m_slices(slices)
  , m_twistPoints(twistPoints)
{
  m_rod.add_property(m_frame, "render frame");
  m_rod.add_property(m_points, "render points");
  m_rod.add_property(m_uvCoord, "render uv-coordinates");
  m_rod.property_handle(m_twist, "twist");
  setSlices(slices);
}

void RodTube::setSlices(int slices)
{
  m_slices = slices;
  computeAngles();
}

int RodTube::getSlices() const
{
  return m_slices;
}

void RodTube::setTwistPoints(bool twistPoints)
{
  m_twistPoints = twistPoints;
}

bool RodTube::getTwistPoints() const
{
  return m_twistPoints;
}

void RodTube::computeAngles()
{
  m_sin2.resize(m_slices);
  m_cos2.resize(m_slices);
  for (int k = 0; k < m_slices; ++k) {
    Scalar angle = 2.0 * M_PI * k / m_slices;
    m_sin2[k] = square(sin(angle));
    m_cos2[k] = square(cos(angle));
  }
}

void RodTube::buildFrame()
{
  ElasticRod::edge_iter eit = m_rod.edges_begin(), end = m_rod.edges_end();
  Vec3d tangent = m_rod.getTangent(*eit);
  Vec3d frame = m_rod.getMaterial1(*eit);
  m_rod.property(m_frame)[*eit] = frame;

  for (++eit; eit != end; ++eit) {
    if (m_twistPoints) {
      m_rod.property(m_frame)[*eit] = m_rod.getMaterial1(*eit);
    } else {
      frame = parallelTransport(frame, tangent, m_rod.getTangent(*eit));
      tangent = m_rod.getTangent(*eit);
      m_rod.property(m_frame)[*eit] = frame;
    }
  }
}

void RodTube::buildTube()
{
  buildFrame();

  Scalar a, b, half, m;
  Vec3d t0, t1, binorm, v, axis;

  Scalar uCoord = 0.0;
  Scalar vCoord = 0.0;
  Scalar scale = m_rod.getRadiusScale() * 1.0;

  for (int i = 0; i < m_rod.nv(); ++i) {
    const Vec3d& vertex = m_rod.getVertex(i);

    if (i == 0) {
      a = scale * m_rod.radiusA(0);
      b = scale * m_rod.radiusB(0);
      axis = m_rod.getTangent(0);
      m = 0;
      v = m_rod.property(m_frame)[0];

    } else if (i == m_rod.nv()-1) {
      a = scale * m_rod.radiusA(i-1);
      b = scale * m_rod.radiusB(i-1);
      axis = m_rod.getTangent(i-1);
      m = 0;
      v = m_rod.property(m_frame)[i-1];

    } else {
      a = scale * 0.5 * (m_rod.radiusA(i-1) + m_rod.radiusA(i));
      b = scale * 0.5 * (m_rod.radiusB(i-1) + m_rod.radiusB(i));

      t0 = m_rod.getTangent(i-1);
      t1 = m_rod.getTangent(i);
      binorm = t0.cross(t1);
      half = 0.5 * atan2(binorm.norm(), t0.dot(t1));
      v = m_rod.property(m_frame)[i-1];
      if( binorm.norm() != 0.0 ) {
        binorm.normalize();
        rotateAxisAngle(v, binorm.normalized(), half);
      }

      axis = (t0+t1).normalized();
      m = m_rod.property(m_twist)[ElasticRod::vertex_handle(i)];
    }

    if (m_twistPoints) {
      rotateAxisAngle(v, axis, m);
      vCoord = 0.0;
    } else {
      vCoord -= m / (2.0 * M_PI);
    }

    v.normalize();
    Scalar aa = a*a;
    Scalar ab = a*b;
    Scalar bb = b*b;
    m_rod.property(m_points)[i].resize(m_slices);
    m_rod.property(m_uvCoord)[i].resize(m_slices);
    for (int k = 0; k < m_slices; ++k) {
      Scalar r = ab / sqrt(bb * m_cos2[k] + aa * m_sin2[k]);
      m_rod.property(m_points)[i][k] = vertex + r * v;
      m_rod.property(m_uvCoord)[i][k] = std::make_pair(uCoord, vCoord);

      rotateAxisAngle(v, axis, 2.0 * M_PI / ((Scalar) m_slices));
      vCoord += 1.0 / m_slices;
    }
    vCoord -= 1.0 + m / (2.0 * M_PI);
    if (i < m_rod.ne()) uCoord += m_rod.getEdgeLength(i);
  }
}

const std::vector<Vec3d>& RodTube::getTubeAtVert(int vertIdx) const
{
  return m_rod.property(m_points)[vertIdx];
}

const Vec3d& RodTube::getPointAtVert(int vertIdx, int pointIdx) const
{
  return m_rod.property(m_points)[vertIdx][pointIdx % m_slices];
}

std::pair<Scalar,Scalar> RodTube::getUVforPointAtVert(int vertIdx, int pointIdx)
{
  int wrap = pointIdx / m_slices;
  std::pair<Scalar, Scalar>& uv
    = m_rod.property(m_uvCoord)[vertIdx][pointIdx % m_slices];
  return std::pair<Scalar, Scalar>(uv.first, uv.second + wrap);
}

} // namespace BASim

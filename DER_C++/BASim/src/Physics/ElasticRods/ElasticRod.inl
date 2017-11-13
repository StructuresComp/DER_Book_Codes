/**
 * \file ElasticRod.inl
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

inline int ElasticRod::edgeIdx(const edge_handle& eh) const
{
  return property(m_edgeIdx)[eh];
}

inline int ElasticRod::edgeIdx(int edgeNumber) const
{
  return property(m_edgeIdx)[edgeNumber];
}

inline const Vec3d& ElasticRod::getVertex(const vertex_handle& vh) const
{
  return property(m_vertexPositions)[vh];
}

inline void ElasticRod::setVertex(const vertex_handle& vh, const Vec3d& v)
{
  property(m_vertexPositions)[vh] = v;
}

inline const Vec3d& ElasticRod::getFromVertex(const edge_handle& eh) const
{
  return getVertex(getEdgeTopology(eh).getFromVertex());
}

inline const Vec3d& ElasticRod::getToVertex(const edge_handle& eh) const
{
  return getVertex(getEdgeTopology(eh).getToVertex());
}

inline const Vec3d& ElasticRod::getVertex(int i) const
{
  return property(m_vertexPositions)[i];
}

inline void ElasticRod::setVertex(int i, const Vec3d& v)
{
  property(m_vertexPositions)[i] = v;
}

inline const Scalar& ElasticRod::getTheta(const edge_handle& eh) const
{
  return property(m_theta)[eh];
}

inline void ElasticRod::setTheta(const edge_handle& eh, const Scalar& t)
{
  property(m_theta)[eh] = t;
}

inline const Scalar& ElasticRod::getTheta(int j) const
{
  return property(m_theta)[j];
}

inline void ElasticRod::setTheta(int j, const Scalar& t)
{
  property(m_theta)[j] = t;
}

inline const Vec3d& ElasticRod::getVelocity(int i) const
{
  return property(m_vertexVelocities)[i];
}

inline void ElasticRod::setVelocity(int i, const Vec3d& v)
{
  property(m_vertexVelocities)[i] = v;
}

inline const Scalar& ElasticRod::getThetaDot(int j) const
{
  return property(m_thetaDot)[j];
}

inline void ElasticRod::setThetaDot(int j, const Scalar& td)
{
  property(m_thetaDot)[j] = td;
}

inline const Vec3d& ElasticRod::getEdge(const edge_handle& eh) const
{
  return property(m_edges)[eh];
}

inline void ElasticRod::setEdge(const edge_handle& eh, const Vec3d& edge)
{
  property(m_edges)[eh] = edge;
}

inline const Vec3d& ElasticRod::getEdge(int j) const
{
  return property(m_edges)[j];
}

inline void ElasticRod::setEdge(int j, const Vec3d& edge)
{
  property(m_edges)[j] = edge;
}

inline const Vec3d& ElasticRod::getTangent(const edge_handle& eh) const
{
  return property(m_tangents)[eh];
}

inline void ElasticRod::setTangent(const edge_handle& eh, const Vec3d& tangent)
{
  property(m_tangents)[eh] = tangent;
}

inline const Vec3d& ElasticRod::getTangent(int j) const
{
  return property(m_tangents)[j];
}

inline void ElasticRod::setTangent(int j, const Vec3d& tangent)
{
  property(m_tangents)[j] = tangent;
}

inline const Scalar& ElasticRod::getVertexMass(int i) const
{
  return property(m_vertexMasses)[i];
}

inline void ElasticRod::setVertexMass(int i, const Scalar& m)
{
  property(m_vertexMasses)[i] = m;
}

inline const Scalar& ElasticRod::getEdgeInertia(int i) const
{
  return property(m_edgeInertias)[i];
}

inline void ElasticRod::setEdgeInertia(int i, const Scalar& I)
{
  property(m_edgeInertias)[i] = I;
}

inline const Vec3d&
ElasticRod::getReferenceDirector1(const edge_handle& eh) const
{
  return property(m_referenceDirectors)[eh].first;
}

inline void
ElasticRod::setReferenceDirector1(const edge_handle& eh, const Vec3d& u)
{
  property(m_referenceDirectors)[eh].first = u;
}

inline const
Vec3d& ElasticRod::getReferenceDirector2(const edge_handle& eh) const
{
  return property(m_referenceDirectors)[eh].second;
}

inline void
ElasticRod::setReferenceDirector2(const edge_handle& eh, const Vec3d& v)
{
  property(m_referenceDirectors)[eh].second = v;
}

inline const Vec3d& ElasticRod::getReferenceDirector1(int i) const
{
  return property(m_referenceDirectors)[i].first;
}

inline void ElasticRod::setReferenceDirector1(int i, const Vec3d& u)
{
  property(m_referenceDirectors)[i].first = u;
}

inline const Vec3d& ElasticRod::getReferenceDirector2(int i) const
{
  return property(m_referenceDirectors)[i].second;
}

inline void ElasticRod::setReferenceDirector2(int i, const Vec3d& v)
{
  property(m_referenceDirectors)[i].second = v;
}

inline const Vec3d& ElasticRod::getMaterial1(const edge_handle& eh) const
{
  return property(m_materialDirectors)[eh].first;
}

inline const Vec3d& ElasticRod::getMaterial2(const edge_handle& eh) const
{
  return property(m_materialDirectors)[eh].second;
}

inline const Vec3d& ElasticRod::getMaterial1(int j) const
{
  return property(m_materialDirectors)[j].first;
}

inline const Vec3d& ElasticRod::getMaterial2(int j) const
{
  return property(m_materialDirectors)[j].second;
}

inline const Scalar& ElasticRod::getEdgeLength(const edge_handle& eh) const
{
  return property(m_edgeLengths)[eh];
}

inline void ElasticRod::setEdgeLength(const edge_handle& eh, const Scalar& len)
{
  property(m_edgeLengths)[eh] = len;
}

inline const Scalar& ElasticRod::getEdgeLength(int j) const
{
  return property(m_edgeLengths)[j];
}

inline void ElasticRod::setEdgeLength(int j, const Scalar& len)
{
  property(m_edgeLengths)[j] = len;
}

inline const Scalar& ElasticRod::getVoronoiLength(int i) const
{
  return property(m_voronoiLengths)[i];
}

inline void ElasticRod::setVoronoiLength(int i, const Scalar& len)
{
  property(m_voronoiLengths)[i] = len;
}

inline bool ElasticRod::quasistatic() const
{
  return property(m_quasistatic);
}

inline void ElasticRod::setQuasistatic(bool q)
{
  property(m_quasistatic) = q;
}

inline ElasticRod::RefFrameType ElasticRod::refFrameType() const
{
  return property(m_refFrameType);
}

inline void ElasticRod::setRefFrameType(RefFrameType type)
{
  property(m_refFrameType) = type;
}

inline const ElasticRod::RodForces& ElasticRod::getForces() const
{
  return property(m_forces);
}

inline ElasticRod::RodForces& ElasticRod::getForces()
{
  return property(m_forces);
}

inline const Scalar& ElasticRod::getYoungsModulus() const
{
  return property(m_YoungsModulus);
}

inline void ElasticRod::setYoungsModulus(const Scalar& E)
{
  property(m_YoungsModulus) = E;
}

inline const Scalar& ElasticRod::getShearModulus() const
{
  return property(m_ShearModulus);
}

inline void ElasticRod::setShearModulus(const Scalar& G)
{
  property(m_ShearModulus) = G;
}

inline const Scalar& ElasticRod::getViscosity() const
{
  return property(m_viscosity);
}

inline void ElasticRod::setViscosity(const Scalar& mu)
{
  property(m_viscosity) = mu;
}

inline int ElasticRod::vertIdx(int vertexNumber, int coordinate) const
{
  return property(m_vertIdx)[vertexNumber] + coordinate;
}

inline int ElasticRod::vertIdx(const vertex_handle& vh, int coordinate) const
{
  return property(m_vertIdx)[vh] + coordinate;
}

inline const Scalar&
ElasticRod::getReferenceTwist(const vertex_handle& vh) const
{
  return property(m_referenceTwist)[vh];
}

inline void
ElasticRod::setReferenceTwist(const vertex_handle& vh,
                              const Scalar& referenceTwist)
{
  property(m_referenceTwist)[vh] = referenceTwist;
}

inline const Scalar& ElasticRod::getReferenceTwist(int i) const
{
  return property(m_referenceTwist)[i];
}

inline void ElasticRod::setReferenceTwist(int i, const Scalar& referenceTwist)
{
  property(m_referenceTwist)[i] = referenceTwist;
}

inline const Vec3d&
ElasticRod::getCurvatureBinormal(const vertex_handle& vh) const
{
  return property(m_curvatureBinormal)[vh];
}

inline void
ElasticRod::setCurvatureBinormal(const vertex_handle& vh, const Vec3d& kb)
{
  property(m_curvatureBinormal)[vh] = kb;
}

inline const Vec3d& ElasticRod::getCurvatureBinormal(int i) const
{
  return property(m_curvatureBinormal)[i];
}

inline void ElasticRod::setCurvatureBinormal(int i, const Vec3d& kb)
{
  property(m_curvatureBinormal)[i] = kb;
}

inline ElasticRod::edge_handle
ElasticRod::inEdge(const vertex_handle& vh) const
{
  const vertex_topology& vt = getVertexTopology(vh);
  return vt[0];
}

inline ElasticRod::edge_handle
ElasticRod::outEdge(const vertex_handle& vh) const
{
  const vertex_topology& vt = getVertexTopology(vh);
  return vt[1];
}

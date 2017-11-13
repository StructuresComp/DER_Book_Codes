/**
 * \file ElasticRod.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef ELASTICROD_HH
#define ELASTICROD_HH

#include "BASim/src/Physics/PhysObject.hh"
#include "BASim/src/Physics/ElasticRods/RodBoundaryCondition.hh"

namespace BASim {

class RodForce;
class MatrixBase;

/** Base class for rods. The degrees of freedom for rods are the
    vertex positions (3 dofs per vertex) and the angles between the
    reference frames and the material frames (1 dof per edge). The
    default indexing that is set up is to interleave the verrtex
    degrees of freedom with the edge degrees of freedom as follows:
    \f$\left<x_0,y_0,z_0,\theta_0,x_1,y_1,z_1,\theta_1,...\right>\f$.
*/
class ElasticRod : public PhysObject
{
public:

  typedef std::vector<RodForce*> RodForces;

  ElasticRod(int numVertices = 3, bool closed = false);
  virtual ~ElasticRod();

	void addVertexAtEnd(Vec3d& x, Vec3d& v);
	void removeVertexFromFront();
	
  /** \name Inherited from PhysObject */

  //@{

  virtual void setup();
  virtual void computeForces(VecXd& force);
  virtual void computeJacobian(int baseidx, Scalar scale, MatrixBase& J);

  virtual const Scalar&
  getVertexDof(const vertex_handle& vh, int num) const;

  virtual void
  setVertexDof(const vertex_handle& vh, int num, const Scalar& dof);

  virtual const Scalar&
  getEdgeDof(const edge_handle& eh, int num) const;

  virtual void
  setEdgeDof(const edge_handle& eh, int num, const Scalar& dof);

  virtual const Scalar&
  getVertexVel(const vertex_handle& vh, int num) const;

  virtual void
  setVertexVel(const vertex_handle& vh, int num, const Scalar& vel);

  virtual const Scalar&
  getEdgeVel(const edge_handle& eh, int num) const;

  virtual void
  setEdgeVel(const edge_handle& eh, int num, const Scalar& vel);

  virtual const Scalar&
  getVertexMass(const vertex_handle& vh, int num) const;

  virtual void
  setVertexMass(const vertex_handle& vh, int num, const Scalar& mass);

  virtual const Scalar&
  getEdgeMass(const edge_handle& eh, int num) const;

  virtual void
  setEdgeMass(const edge_handle& eh, int num, const Scalar& mass);

  //@}

  edge_handle inEdge(const vertex_handle& vh) const;
  edge_handle outEdge(const vertex_handle& vh) const;

  const Scalar& getYoungsModulus() const;
  void setYoungsModulus(const Scalar& E);

  const Scalar& getShearModulus() const;
  void setShearModulus(const Scalar& G);

  const Scalar& getViscosity() const;
  void setViscosity(const Scalar& mu);

  virtual int vertIdx(const vertex_handle& vh, int coordinate) const;
  virtual int edgeIdx(const edge_handle& eh) const;

  virtual int vertIdx(int vertexNumber, int coordinate) const;
  virtual int edgeIdx(int edgeNumber) const;

  const Vec3d& getVertex(const vertex_handle& vh) const;
  void setVertex(const vertex_handle& vh, const Vec3d& v);

  const Vec3d& getFromVertex(const edge_handle& eh) const;
  const Vec3d& getToVertex(const edge_handle& eh) const;

  const Vec3d& getVertex(int i) const;
  void setVertex(int i, const Vec3d& v);

  const Scalar& getTheta(const edge_handle& eh) const;
  void setTheta(const edge_handle& eh, const Scalar& t);

  const Scalar& getTheta(int j) const;
  void setTheta(int j, const Scalar& t);

  const Vec3d& getVelocity(int i) const;
  void setVelocity(int i, const Vec3d& v);

  const Scalar& getThetaDot(int j) const;
  void setThetaDot(int j, const Scalar& td);

  const Vec3d& getEdge(const edge_handle& eh) const;
  void setEdge(const edge_handle& eh, const Vec3d& edge);

  const Vec3d& getEdge(int j) const;
  void setEdge(int j, const Vec3d& edge);

  const Vec3d& getTangent(const edge_handle& eh) const;
  void setTangent(const edge_handle& eh, const Vec3d& tangent);

  const Vec3d& getTangent(int j) const;
  void setTangent(int j, const Vec3d& tangent);

  const Scalar& getVertexMass(int i) const;
  void setVertexMass(int i, const Scalar& m);

  const Scalar& getEdgeInertia(int j) const;
  void setEdgeInertia(int j, const Scalar& I);

  /** Sets a circular cross-section for all of the edges in the rod.

      \param[in] r The radius of the cross-section.
  */
  void setRadius(const Scalar& r);

  /** Sets an elliptical cross-section for all of the edges in the rod.

      \param[in] a The major axis of the ellipse.
      \param[in] b The minor axis of the ellipse.
  */
  void setRadius(const Scalar& a, const Scalar& b);

  /** Sets a circular cross-section for a specific edge in the rod.

      \param[in] j The edge number.
      \param[in] r The radius of the cross-section.
  */
  void setRadius(int j, const Scalar& r);
  void setRadius(const edge_handle& eh, const Scalar& r);

  /** Sets an elliptical cross-section for a specific edge in the rod.

      \param[in] j The edge number.
      \param[in] a The major radius of the ellipse.
      \param[in] b The minor radius of the ellipse.
  */
  void setRadius(int j, const Scalar& a, const Scalar& b);

  /** Gets the radius of the rod. It only makes sense to use this call
      if all of the edges in the rod have uniform, circular
      cross-sections.

      \return The radius of the rod.
  */
  const Scalar& radius() const;

  /** Gets the radius of an edge in the rod. It only makes sense to
      use this call if the edge has a circular cross-section.

      \param[in] j The edge number.
      \return The radius of the edge.
  */
  const Scalar& radius(int j) const;

  /** Gets the major radius of an edge in the rod.

      \param[in] j The edge number.
      \return The major radius of the edge.
  */
  const Scalar& radiusA(int j) const;
  const Scalar& radiusA(const edge_handle& eh) const;

  /** Gets the minor radius of an edge in the rod.

      \param[in] j The edge number.
      \return The minor radius of the edge.
  */
  const Scalar& radiusB(int j) const;
  const Scalar& radiusB(const edge_handle& eh) const;

  /** Multiplicative factor that scales the radius for rendering and
      collision detection. Defaults to 1.

      \return radius scale factor
   */
  Scalar getRadiusScale() const;
  void setRadiusScale(Scalar s);

  const Vec3d& getReferenceDirector1(const edge_handle& eh) const;
  void setReferenceDirector1(const edge_handle& eh, const Vec3d& u);
  const Vec3d& getReferenceDirector2(const edge_handle& eh) const;
  void setReferenceDirector2(const edge_handle& eh, const Vec3d& v);

  const Vec3d& getReferenceDirector1(int j) const;
  void setReferenceDirector1(int j, const Vec3d& u);
  const Vec3d& getReferenceDirector2(int j) const;
  void setReferenceDirector2(int j, const Vec3d& v);

  const Vec3d& getMaterial1(const edge_handle& eh) const;
  void setMaterial1(const edge_handle& eh, const Vec3d& m1);
  const Vec3d& getMaterial2(const edge_handle& eh) const;
  void setMaterial2(const edge_handle& eh, const Vec3d& m2);

  const Vec3d& getMaterial1(int j) const;
  void setMaterial1(int j, const Vec3d& m1);
  const Vec3d& getMaterial2(int j) const;
  void setMaterial2(int j, const Vec3d& m2);

  const Scalar& getEdgeLength(const edge_handle& eh) const;
  void setEdgeLength(const edge_handle& eh, const Scalar& len);

  const Scalar& getEdgeLength(int j) const;
  void setEdgeLength(int j, const Scalar& len);

  const Scalar& getVoronoiLength(int i) const;
  void setVoronoiLength(int i, const Scalar& len);

  const Scalar& getReferenceTwist(const vertex_handle& vh) const;
  void setReferenceTwist(const vertex_handle& vh, const Scalar& referenceTwist);

  const Scalar& getReferenceTwist(int i) const;
  void setReferenceTwist(int i, const Scalar& referenceTwist);

  const Vec3d& getCurvatureBinormal(const vertex_handle& vh) const;
  void setCurvatureBinormal(const vertex_handle& vh, const Vec3d& kb);

  const Vec3d& getCurvatureBinormal(int i) const;
  void setCurvatureBinormal(int i, const Vec3d& kb);

  bool quasistatic() const;
  void setQuasistatic(bool q);

  enum RefFrameType { SpaceParallel, TimeParallel };
  RefFrameType refFrameType() const;
  void setRefFrameType(RefFrameType type);

  const Scalar& density() const { return property(m_density); }
  void setDensity(const Scalar& d) { property(m_density) = d; }

  void setHollow(const bool& h, const Scalar& r) {
    property(m_is_hollow) = h;
    property(m_inner_radius) = r;
  }
  
  const bool& getIsHollow() const { return property(m_is_hollow); }
  const Scalar& getInnerRadius() const { return property(m_inner_radius); }

  const RodForces& getForces() const;
  RodForces& getForces();
  void addForce(RodForce* force);

  RodBoundaryCondition* getBoundaryCondition()
  {
    if (m_boundaryConditions == NULL) m_boundaryConditions = new RodBoundaryCondition(*this);
    return m_boundaryConditions;
  }

  const RodBoundaryCondition* getBoundaryCondition() const
  {
    return m_boundaryConditions;
  }

  virtual void updateProperties();
  virtual void updateReferenceProperties();
  virtual void verifyProperties();

  /** At the beginning of the time step, the undeformed configuration
      must be set to the current configuration for the viscous
      forces */
  virtual void viscousUpdate();

  /** The stiffness of the viscous forces (internal damping) depend on
      the size of the time step being taken, so they must be
      recomputed whenever the size of the time step is changed. */
  void setTimeStep(Scalar dt);
  Scalar getTimeStep() const { return property(m_dt); }


  void computeEdges();
  void computeTangents();
  void computeEdgeLengths();
  void computeVoronoiLengths();
  void computeReferenceDirectors();
  void computeSpaceParallel();
  void computeTimeParallel();
  void computeCurvatureBinormals();
  void computeReferenceTwist();
  void computeMaterialDirectors();
  void computeVertexMasses();
  void computeEdgeInertias();
  void setupDofIndices();

  void updateForceProperties();

  VecXd m_twisting_energy;

  void setBendingStiffnessMultiplier(int v, Scalar m) { property(m_bendingStiffnessMultiplier)[v] = m; }
  Scalar getBendingStiffnessMultiplier(int v) { return property(m_bendingStiffnessMultiplier)[v]; }
  void setBendingStiffnessMultiplier(const vertex_handle & v, Scalar m) { property(m_bendingStiffnessMultiplier)[v] = m; }
  Scalar getBendingStiffnessMultiplier(const vertex_handle & v) { return property(m_bendingStiffnessMultiplier)[v]; }
  
  void setViscousBendingStiffnessOverride(int v, Scalar m) { property(m_viscousBendingStiffnessOverride)[v] = m; }
  Scalar getViscousBendingStiffnessOverride(int v) { return property(m_viscousBendingStiffnessOverride)[v]; }
  void setViscousBendingStiffnessOverride(const vertex_handle & v, Scalar m) { property(m_viscousBendingStiffnessOverride)[v] = m; }
  Scalar getViscousBendingStiffnessOverride(const vertex_handle & v) { return property(m_viscousBendingStiffnessOverride)[v]; }
  
protected:

  /** Computes the mass of an elliptical cylinder, which is the
      generic representation of each edge of the rod.

      \param[in] density The volumetric density of the cylinder.
      \param[in] a The major radius of the cylinder.
      \param[in] b The minor radius of the cylinder.
      \param[in] h The height of the cylinder.
      \return The mass of the cylinder.
  */
  Scalar computeMass(Scalar density, Scalar a, Scalar b, Scalar h);

  ObjPropHandle<RodForces> m_forces; ///< forces acting on the rod
  ObjPropHandle<bool> m_quasistatic;
  ObjPropHandle<RefFrameType> m_refFrameType;
  ObjPropHandle<Scalar> m_density;
  ObjPropHandle<Scalar> m_YoungsModulus;
  ObjPropHandle<Scalar> m_ShearModulus;
  ObjPropHandle<Scalar> m_viscosity;
  ObjPropHandle<Scalar> m_dt;
  ObjPropHandle<Scalar> m_radius_scale;

  ObjPropHandle<bool> m_is_hollow;
  ObjPropHandle<Scalar> m_inner_radius;

  VPropHandle<Vec3d> m_vertexPositions;
  VPropHandle<Vec3d> m_vertexVelocities;
  VPropHandle<Scalar> m_voronoiLengths;
  VPropHandle<Scalar> m_vertexMasses;
  VPropHandle<Scalar> m_referenceTwist; ///< twist of the reference frame
  VPropHandle<Vec3d> m_curvatureBinormal;
  VPropHandle<int> m_vertIdx;

  EPropHandle<Scalar> m_theta;
  EPropHandle<Scalar> m_thetaDot;
  EPropHandle< Util::pair<Scalar, Scalar> > m_edgeRadius;
  EPropHandle<Scalar> m_edgeInertias;
  EPropHandle< Util::pair<Vec3d,Vec3d> > m_referenceDirectors;
  EPropHandle< Util::pair<Vec3d,Vec3d> > m_materialDirectors;
  EPropHandle<Vec3d> m_edges;
  EPropHandle<Vec3d> m_tangents;
  EPropHandle<Scalar> m_edgeLengths; ///< lengths of edges
  EPropHandle<int> m_edgeIdx;
  
  VPropHandle<Scalar> m_bendingStiffnessMultiplier;
  VPropHandle<Scalar> m_viscousBendingStiffnessOverride;

  RodBoundaryCondition* m_boundaryConditions;
};

#include "ElasticRod.inl"

} // namespace BASim

#endif // ELASTICROD_HH

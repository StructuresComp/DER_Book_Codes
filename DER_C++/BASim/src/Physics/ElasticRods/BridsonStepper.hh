/**
 * \file BridsonStepper.hh
 *
 * \author smith@cs.columbia.edu
 * \date 02/16/2010
 */

// This class is currently a little hacked together, and
// does some inefficient stuff (copies all DOFs 
// into flat arrays and such) just because it was
// easier than using BASim iterators, for now :).

#ifndef BRIDSONSTEPPER_HH
#define BRIDSONSTEPPER_HH

#include "BASim/src/Core/ObjectControllerBase.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Physics/ElasticRods/RodTimeStepper.hh"
#include "BASim/src/Physics/ElasticRods/RodMassDamping.hh"
#include "BASim/src/Physics/ElasticRods/RodGravity.hh"

#include "BASim/src/Physics/ElasticRods/RodDampingFriction.hh"

#include "BASim/src/Math/Math.hh"

#include "BASim/src/Core/TriangleMesh.hh"

#include "BASim/src/Collisions/BVHAABB.hh"
#include "BASim/src/Collisions/CollisionUtils.hh"

#include "MinimalRodStateBackup.hh"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <limits>

namespace BASim {

  /**
   * Class to evolve a collection of rods forward in time, resolving collisions using
   * a "velocity filter" in the spirit of Bridson's 2002 paper "Robust Treatment of
   * Collisions, Contact, and Friction for Cloth Animation."
   */
  class BridsonStepper : public ObjectControllerBase
  {

  public:
    /**
     * Default constructor.
     */
    BridsonStepper();

    /**
     * Creates a BridsonStepper with user-supplied options. 
     *
     * \param[in] intgrtr Integrator (class RodTimeStepper) to use. Assumes implicit euler, for now.
     * \param[in] max_implct_itrtns If an implicit integrator is selected, the maximum iterations allowed.
     * \param[in] dt Timestep to use.
     * \param[in] mass_dmpng Amount of damping that acts in opposition to vertex velocities (I think? Miklos has been mucking with the damping :)).
     * \param[in] grav Three dimensional vector that specifies gravity.
     */
    BridsonStepper( const RodTimeStepper::Method& intgrtr, const int& max_implct_itrtns, const double& dt, const double& mass_dmpng, const Vec3d& grav );

    /**
     * Destructor.
     */
    virtual ~BridsonStepper();
    
    /**
     * After adding new rods or objects, this method must be called.
     */
    void prepareForExecution();
    
    /**
     * Adds a rod that will be evolved in time using this BridsonStepper.
     */
    void addRod( ElasticRod* rod );
    
    /**
     * Adds a non-simulated triangle mesh for objects to collide with.
     */
    void addTriangleMesh( TriangleMesh* tri_mesh );

    /**
     * Evolves all inserted rods forward in time. 
     */
    void execute();

    void execute2();
    
    /**
     * Modifies the timestep.
     */
    void setDt( const double& dt )
    {
      assert( dt > 0.0 );
      m_dt = dt;
    }

    /**
     * An unneccessary (sp?) method that I am deleting soon :). 
     */
    void setGravity( const Vec3d& gravity ) { m_gravity = gravity; };

    /**
     * Disables all response.
     */
    void disableResponse();

    /**
     * Enables response, subject to state of enablePenaltyImpulses() and enableIterativeInelasticImpulses().
     */
    void enableReseponse();

    /**
     * Enables penalty response.
     */
    void enablePenaltyImpulses();

    /**
     * Disables penalty response.
     */
    void disablePenaltyImpulses();

    /**
     * Enables iterative impulse response.
     */
    void enableIterativeInelasticImpulses();

    /**
     * Disables iterative impulse response.
     */
    void disableIterativeInelasticImpulses();

    /**
     * Sets the maximum number of inelastic impulses to apply iterativly.
     */
    void setNumInelasticIterations( const int& num_itr );

    /**
     * Sets the stiffness of edge-edge penalty impulses. 
     */
    void setEdgeEdgePenalty( const double& k );

    /**
     * Sets the stiffness of vertex-face penalty impulses. 
     */
    void setVertexFacePenalty( const double& k );

    /**
     * Number of rods this controller is responsible for.
     */
    int getNumRods() const { return m_rods.size(); };

    /**
     * Number of triangle meshes this controller is responsible for.
     */
    int getNumTriangleMeshes() const { return m_triangle_meshes.size(); };
    
    // TODO: Move these to some kind of automated test suite    
    void testCoplanarityTime();
    
    RodTimeStepper* getInternalRodTimeStepper(int i);
    
    void resizeRodVertex();
    
    void setContactAndFrictionMethod( int contact, int friction ) 
    {
			m_contact_method = contact;
			m_friction_method = friction;    	
    }
    
    void setFriction( bool f );

    void setCoulombFriction( double static_cof, double kinetic_cof );
    
    void setGroundFix( bool flag ) {
    	m_ground_fix = flag;
    }

    void setTwistFixAfterHit( bool flag ) {
    	m_fix_material_frame = flag;
    }
    
    void removeRod();
    
    
    
  private:
    /////////////////////////////////////////////////////
    // Methods for checking the sanity of input rods

    // Currently we do not support collisions for anisotropic cross sections
    void ensureCircularCrossSection( const ElasticRod& rod ) const;

    // If the cross-sectional radius is too large and edge lengths too small,
    // non-adjacent portions of the rod will be in contact by default. We can
    // probably add some special case code to handle this later, but just 
    // disallow the situation for now.
    void ensureNoCollisionsByDefault( const ElasticRod& rod ) const;


    /////////////////////////////////////////////////////
    // Helper methods

    // Returns the total number of degrees of freedom in the system
    int getNumDof() const;
    
    // Returns the total number of vertices in the system
    int getNumVerts() const;
    
    // Creates a time stepper for a given rod
    RodTimeStepper* createRodTimeStepper(ElasticRod& rod);

    void extractPositions( const std::vector<ElasticRod*>& rods, const std::vector<int>& base_indices, VecXd& positions );
    void extractVelocities( const std::vector<ElasticRod*>& rods, const std::vector<int>& base_indices, VecXd& velocities );
    //void extractMasses( const std::vector<ElasticRod*>& rods, std::vector<double>& masses ) const;
    //void extractFixedVertices( const std::vector<ElasticRod*>& rods, std::vector<bool>& fixed ) const;
    //void extractEdges( const std::vector<ElasticRod*>& rods, const std::vector<int>& base_indices, std::vector<std::pair<int,int> >& edges, std::vector<double>& radii );

//		void extractTheta( const std::vector<ElasticRod*>& rods, VecXd& angles );
		void extractTheta( const std::vector<ElasticRod*>& rods, VecXd& angles, VecXd& angleDots );
    
    void restorePositions( std::vector<ElasticRod*>& rods, const VecXd& positions );
    void restoreVelocities( std::vector<ElasticRod*>& rods, const VecXd& velocities );
		void restoreTheta( const std::vector<ElasticRod*>& rods, const VecXd& angles, const VecXd& angleDots );
    

    /////////////////////////////////////////////////////
    // Collision detection routines
    
    // Determines if two edges share a vertex
    bool edgesShareVertex( const std::pair<int,int>& edgei, const std::pair<int,int>& edgej ) const;
    bool edgesSharevertex( const int& e0v0, const int& e0v1, const int& e1v0, const int& e1v1 ) const;

    
    // Determines if a vertex and a face share a vertex
    bool vertexAndFaceShareVertex( const int& vertex, const int& face ) const;
    bool vertexAndFaceShareVertex( const int& v, const int& f0, const int& f1, const int& f2 ) const;

    
    // Generates a list of ALL possible edge-edge collisions
    void generateAllEdgeEdgeProximityCollisionPairs( std::vector<EdgeEdgeProximityCollision>& edge_edge_collisions ) const;
    void generateAllEdgeEdgeContinuousTimeCollisionPairs( std::vector<EdgeEdgeContinuousTimeCollision>& edge_edge_collisions ) const;

    // Generates a list of ALL possible vertex-face collisions
    void generateAllVertexFaceProximityCollisionPairs( std::vector<VertexFaceProximityCollision>& vertex_face_collisions ) const;
    void generateAllVertexFaceContinuousTimeCollisionPairs( std::vector<VertexFaceContinuousTimeCollision>& vertex_face_collisions ) const;

    // Returns a list of all edges that are in "proximity" for the positions specified in the vector x.
    //void detectEdgeEdgeProximityCollisions( const VecXd& x, std::vector<EdgeEdgeProximityCollision>& edge_edge_collisions ) const;
    void detectEdgeEdgeProximityCollisions( const VecXd& x, std::vector<EdgeEdgeProximityCollision>& pssbl_cllsns, std::vector<EdgeEdgeProximityCollision>& cllsns ) const;

    // Returns a list of all vertex-face pairs that are in "proximity" for the positions specified in the vector x.
    void detectVertexFaceProximityCollisions( const VecXd& x, std::vector<VertexFaceProximityCollision>& pssbl_cllsns, std::vector<VertexFaceProximityCollision>& vetex_face_collisions ) const;

    // Detects all edge-edge continuous time collisions
    void detectEdgeEdgeContinuousTimeCollisions( const VecXd& x, const VecXd& v, std::vector<EdgeEdgeContinuousTimeCollision>& edge_edge_collisions, std::vector<EdgeEdgeContinuousTimeCollision>& edge_edge_collisions_for_real );

    // Detects all vertex-face continuous time collisions
    void detectVertexFaceContinuousTimeCollisions( const VecXd& x, const VecXd& v, std::vector<VertexFaceContinuousTimeCollision>& vertex_face_collisions, std::vector<VertexFaceContinuousTimeCollision>& vertex_face_collisions_for_real );

    // Computes the relative velocity of fixed pieces of material (s and t assumed constant) of two edges
    Vec3d computeRelativeVelocity( const VecXd& v, const int& idxa0, const int& idxa1, const int& idxb0, const int& idxb1, const double& s, const double& t );

    // Computes the relative velocity of fixed pieces of material (u, v, and w assumed constant) of a vertex and face
    Vec3d computeRelativeVelocity( const VecXd& vel, const int& vrtidx, const int& fcidx0, const int& fcidx1, const int& fcidx2, const double& u, const double& v, const double& w );

    // Computes the impulse necessary to eliminate all relative velocity at given points on two edges
    Vec3d computeEdgeEdgeInelasticImpulse( const double& ma0, const double& ma1, const double& mb0, const double& mb1,
                                           const double& s, const double& t, const double& relvel, const Vec3d& n );
    
    // Computes the impulse necessary to eliminate all relative velocity at a given point on a face and a vertex
    Vec3d computeVertexFaceInelasticImpulse( const double& mvrt, const double& mfc0, const double& mfc1, const double& mfc2,
                                             const double& u, const double& v, const double& w, const double& relvel, const Vec3d& n );

    /////////////////////////////////////////////////////
    // Collision response routines
    
    void exertPenaltyImpulses( std::vector<EdgeEdgeProximityCollision>& edg_edg_cllsns, std::vector<VertexFaceProximityCollision>& vrtx_fce_cllsns, VecXd& v );
    void exertInelasticImpulses( std::vector<EdgeEdgeContinuousTimeCollision>& edg_edg_cllsns, std::vector<VertexFaceContinuousTimeCollision>& vrtx_fce_cllsns, VecXd& v );
    
    void exertVertexImpulse( const Vec3d& I, const double& m, const int& idx, VecXd& v );
    void exertEdgeImpulse( const Vec3d& I, const double& m0, const double& m1, const double& alpha, const int& idx0, const int& idx1, VecXd& v );
    void exertFaceImpulse( const Vec3d& I, const double& m0, const double& m1, const double& m2, 
                           const double& u, const double& v, const double& w, 
                           const int& idx0, const int& idx1, const int& idx2, VecXd& vel );
    
    void exertVertexFrictionImpulse( const Vec3d& I, const double& m, const int& idx, VecXd& v, Vec3d& rel_vel, Vec3d& normal);
    
    
    /////////////////////////////////////////////////////////////////
    // Collision detection code adapted from Robert Bridson's website
    
    void addUnique( std::vector<double>& a, double e ) const;
    double triple( const Vec3d &a, const Vec3d &b, const Vec3d &c ) const;
    double signed_volume( const Vec3d& x0, const Vec3d& x1, const Vec3d& x2, const Vec3d& x3 ) const;
    void getCoplanarityTimes( const Vec3d& x0, const Vec3d& x1, const Vec3d& x2, const Vec3d& x3,
                              const Vec3d& xnew0, const Vec3d& xnew1, const Vec3d& xnew2, const Vec3d& xnew3,
                              std::vector<double>& times, std::vector<double>& errors ) const;
    




    // Total number of degrees of freedom in the system
    int m_num_dof;
    // Vector of rods this BridsonStepper evolves in time
    std::vector<ElasticRod*> m_rods;
    // Vector of ScriptedTriangleObjects in the system
    std::vector<TriangleMesh*> m_triangle_meshes;
    // Time steppers to evolve rods forward (ignoring collisions)
    std::vector<RodTimeStepper*> m_steppers;
    // Integrator selected by user
    RodTimeStepper::Method m_method;
    // Timestep selected by user
    double m_dt;
    // Drag set by user
    double m_mass_damping;
    // Gravity selected by user
    Vec3d m_gravity;
    // Number of iterations of implicit solver selected by user
    int m_max_implicit_iterations;
    
    // Entry i is base index of ith rod in global array of position dofs
    std::vector<int> m_base_indices;

    // Entry i is base index of ith ScriptedTriangleMesh in global array of position dofs
    std::vector<int> m_base_triangle_indices;

    // Vector of edges in the system. FREE (not part of a face) edges
    std::vector<std::pair<int,int> > m_edges;
    
    // Vector of triangular faces in the system
    std::vector<TriangularFace> m_faces;

    std::vector<double> m_vertex_radii;
    // TODO: Possibly get rid of these, just pull radii from m_vertex_radii.
    std::vector<double> m_edge_radii;
    std::vector<double> m_face_radii;
    
    std::vector<double> m_masses;
    
    VecXd m_xn;
    VecXd m_xnp1;
    VecXd m_vnphalf;
    
    // Enable/Disable portions of the collision response
    bool m_respns_enbld;
    bool m_pnlty_enbld;
    bool m_itrv_inlstc_enbld;
    int m_num_inlstc_itrns;
    double m_edg_edg_pnlty;
    double m_vrt_fc_pnlty;
    
    bool m_friction_enbld;
    double m_friction_cof;
    
    double m_coulomb_static;
    double m_coulomb_kinetic;
    
    bool m_ground_fix;
    
    // An external force for friction
    RodDampingFriction* m_friction_force;

    // Some debug stuff in for now.
    bool m_nan_enc;
    bool m_inf_enc;
    bool m_lt0_enc;
    bool m_gt0_enc;
    
    BVHAABB* m_bvh;
    
    bool m_fix_material_frame;

		int m_contact_method;
		int m_friction_method;    
		
		MinimalRodStateBackup* rodBackup;

  };
  
} // namespace BASim

#endif // RODTIMESTEPPER_HH


// TODO:
//  o Bridson 2002 does some extra trapezoid-rule stuff at end of step. Look into this.
//  o Add a little bit of extra "kick" to (currently) inelastic impulses to help with FPA errors
//  o Clean up the internal state a bit
//  o Add support for remaining corner cases in detection
//  o There was some subtle indexing bug that was counting collisions twice, add some regression tests for this :). 
//    The check of separting velocities ensures we don't overshoot, but a little caution won't hurt.

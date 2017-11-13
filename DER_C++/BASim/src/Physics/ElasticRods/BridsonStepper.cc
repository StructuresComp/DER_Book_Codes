/**
 * \file BridsonStepper.cc
 *
 * \author smith@cs.columbia.edu
 * \date 02/16/2010
 */

#include "BridsonStepper.hh"


namespace BASim 
{
  
BridsonStepper::BridsonStepper()
:m_num_dof(0)
,m_rods()
,m_triangle_meshes()
,m_steppers()
,m_method(RodTimeStepper::IMPL_EULER)
,m_dt(0.01)
,m_mass_damping(0.0)
,m_gravity(0.0,0.0,0.0)
,m_max_implicit_iterations(10000000) 
,m_base_indices()
,m_base_triangle_indices()
,m_edges()
,m_faces()
,m_respns_enbld(true)
,m_pnlty_enbld(true)
,m_itrv_inlstc_enbld(true)
,m_num_inlstc_itrns(10000)
,m_edg_edg_pnlty(200.0)
,m_vrt_fc_pnlty(200.0)
,m_friction_enbld(false)
,m_friction_cof(0.0)
,m_coulomb_static(0.36)
,m_coulomb_kinetic(0.31)
,m_ground_fix(false)
,m_friction_force(NULL)
,m_nan_enc(false)
,m_inf_enc(false)
,m_lt0_enc(false)
,m_gt0_enc(false)
,m_bvh(NULL)
,m_fix_material_frame(false)
,rodBackup(NULL)
{}


BridsonStepper::BridsonStepper( const RodTimeStepper::Method& intgrtr, const int& max_implct_itrtns, const double& dt, const double& mass_dmpng, const Vec3d& grav )
:m_num_dof(0)
,m_rods()
,m_triangle_meshes()
,m_steppers()
,m_method(intgrtr)
,m_dt(dt)
,m_mass_damping(mass_dmpng)
,m_gravity(grav)
,m_max_implicit_iterations(max_implct_itrtns)
,m_base_indices()
,m_base_triangle_indices()
,m_edges()
,m_faces()
,m_respns_enbld(true)
,m_pnlty_enbld(true)
,m_itrv_inlstc_enbld(true)
,m_num_inlstc_itrns(10000)
,m_edg_edg_pnlty(200.0)
,m_vrt_fc_pnlty(200.0)
,m_friction_enbld(false)
,m_friction_cof(0.0)
,m_coulomb_static(0.36)
,m_coulomb_kinetic(0.31)
,m_ground_fix(false)
,m_friction_force(NULL)
,m_nan_enc(false)
,m_inf_enc(false)
,m_lt0_enc(false)
,m_gt0_enc(false)
,m_bvh(NULL)
,m_fix_material_frame(false)
,rodBackup(NULL)
{
  assert( dt > 0.0 );
  assert( mass_dmpng >= 0.0 );
  assert( m_method == RodTimeStepper::IMPL_EULER );
  if( m_method == RodTimeStepper::IMPL_EULER ) assert( max_implct_itrtns > 0 );
  //m_bvh = new BVHAABB( m_xn, m_edges, m_faces, m_vertex_radii );
  
  rodBackup = new MinimalRodStateBackup();
}


BridsonStepper::~BridsonStepper()
{
  assert( m_rods.size() == m_steppers.size() );
  
  for( int i = 0; i < (int) m_steppers.size(); ++i )
  {
    assert( m_steppers[i] != NULL );
    delete m_steppers[i];
    m_steppers[i] = NULL;
  }
  
  if( m_bvh != NULL )
  {
    delete m_bvh;
    m_bvh = NULL;
  }
  
  if (rodBackup)
  { 
  	delete rodBackup;
  	rodBackup = NULL;
  }
}

void BridsonStepper::removeRod()
{
  m_rods.clear();
  resizeRodVertex();
}

void BridsonStepper::prepareForExecution()
{
  if( m_bvh != NULL )
  {
    delete m_bvh;
    m_bvh = NULL;
  }

  if( m_respns_enbld )
  {
    // Load positions for initial construction of the BVH
    extractPositions(m_rods,m_base_indices,m_xn);
    extractVelocities(m_rods,m_base_indices,m_vnphalf);

    m_bvh = new BVHAABB( m_xn, m_vnphalf, m_edges, m_faces, m_vertex_radii, m_dt );
  }
}

void BridsonStepper::addRod( ElasticRod* rod )
{
  assert( rod != NULL );
  
  // Sanity checks for collision detection purpouses
  ensureCircularCrossSection( *rod );
  ensureNoCollisionsByDefault( *rod );

  // Extract edges from the new rod
  for( int i = 0; i < (int) rod->nv()-1; ++i )
  {
    m_edges.push_back(std::pair<int,int>( getNumVerts()+i, getNumVerts()+i+1 ) );        
    assert( rod->getRadiusScale()*rod->radiusA(i) > 0.0 );
    m_edge_radii.push_back(rod->getRadiusScale()*rod->radiusA(i));
  }
  assert( m_edges.size() == m_edge_radii.size() );

  for( int i = 0; i < (int) rod->nv()-1; ++i )
  {
    assert( rod->getRadiusScale()*rod->radiusA(i) > 0.0 );
    // Radii associated with edges ... what to do if at a vertex with two radii? Average?
    m_vertex_radii.push_back(rod->getRadiusScale()*rod->radiusA(i));
  }
  m_vertex_radii.push_back(m_vertex_radii.back());
  
  // Extract masses from the new rod
  for( ElasticRod::vertex_iter itr = rod->vertices_begin(); itr != rod->vertices_end(); ++itr )
  {
    if( rod->getBoundaryCondition()->isVertexScripted((*itr).idx()) ) 
    {
      m_masses.push_back(std::numeric_limits<double>::infinity());
    }
    else
    {
      assert( rod->getVertexMass(*itr,-1) > 0.0 );
      m_masses.push_back( rod->getVertexMass(*itr,-1) );
    }
  }
  
  // Add the new rod to the system
  m_rods.push_back(rod);

  // Create a timestepper to evolve the rod in time
  RodTimeStepper* newstepper = this->createRodTimeStepper(*rod);
  m_steppers.push_back(newstepper);

  // Update vector that tracks the rod DOF in the system
  m_base_indices.push_back(getNumDof());
  
  // Update total number of DOF in the system
  m_num_dof += 3*rod->nv();
  
  // Resize the internal temporary storage
  m_xn.resize(getNumDof()); m_xn.setConstant(std::numeric_limits<double>::signaling_NaN());
  m_xnp1.resize(getNumDof()); m_xnp1.setConstant(std::numeric_limits<double>::signaling_NaN());
  m_vnphalf.resize(getNumDof()); m_vnphalf.setConstant(std::numeric_limits<double>::signaling_NaN());

  //m_bvh->regenerateBVH( m_xn, m_edges, m_faces, m_vertex_radii );
  //m_bvh->regenerateBVH();
  
  assert( m_rods.size() == m_steppers.size() );
  assert( m_rods.size() == m_base_indices.size() );
  assert( (int) m_masses.size() == getNumVerts() );
  assert( (int) m_vertex_radii.size() == getNumVerts() );
  
//  resizeRodVertex();
}

void BridsonStepper::setFriction( bool f ) {

  m_friction_enbld = f;
  
  if (f) {
   	std::cout << "Friction turned on\n";
  
    if (m_friction_force) {
    	delete m_friction_force;
    }
  
    m_friction_force = new RodDampingFriction( );
    m_friction_force->setContactAndFrictionMethod( m_contact_method, m_friction_method);
    m_friction_force->enable();
  
    m_friction_force->setCoulombCoefficient( m_coulomb_static, m_coulomb_kinetic );
  }		

}

void BridsonStepper::setCoulombFriction( double static_cof, double kinetic_cof ) {
  
  m_coulomb_static = static_cof;
  m_coulomb_kinetic = kinetic_cof;
      
  m_friction_force->setCoulombCoefficient( static_cof, kinetic_cof );

}

void BridsonStepper::resizeRodVertex() {

	m_base_indices.clear();
	m_base_triangle_indices.clear();
	m_edges.clear();
	m_faces.clear();
	m_vertex_radii.clear();
	m_edge_radii.clear();
	m_face_radii.clear();
	m_masses.clear();

	m_num_dof = 0;
	
	for(int r=0; r<(int)m_rods.size(); r++) {
		ElasticRod* rod = m_rods[r];
		
		for( int i = 0; i < (int) rod->nv()-1; ++i )
		{
		  m_edges.push_back(std::pair<int,int>( getNumVerts()+i, getNumVerts()+i+1 ) );        
		  assert( rod->getRadiusScale()*rod->radiusA(i) > 0.0 );
		  m_edge_radii.push_back(rod->getRadiusScale()*rod->radiusA(i));
		}	
		assert( m_edges.size() == m_edge_radii.size() );

		for( int i = 0; i < (int) rod->nv()-1; ++i )
		{
		  assert( rod->getRadiusScale()*rod->radiusA(i) > 0.0 );
		  // Radii associated with edges ... what to do if at a vertex with two radii? Average?
		  m_vertex_radii.push_back(rod->getRadiusScale()*rod->radiusA(i));
		}
		m_vertex_radii.push_back(m_vertex_radii.back());	
		
		// Extract masses from the new rod
		for( ElasticRod::vertex_iter itr = rod->vertices_begin(); itr != rod->vertices_end(); ++itr )
		{
		  if( rod->getBoundaryCondition()->isVertexScripted((*itr).idx()) ) 
		  {
		    m_masses.push_back(std::numeric_limits<double>::infinity());
		  }
		  else
		  {
		    assert( rod->getVertexMass(*itr,-1) > 0.0 );
		    m_masses.push_back( rod->getVertexMass(*itr,-1) );
		  }
		}		

		// Update vector that tracks the rod DOF in the system
		m_base_indices.push_back(getNumDof());
		
		// Update total number of DOF in the system
		m_num_dof += 3*rod->nv();	
	}	


	for(int r=0; r<(int)m_triangle_meshes.size(); r++) {
		TriangleMesh* tri_mesh = m_triangle_meshes[r];
		
		bool fixed = tri_mesh->fixed;

		// Extract faces from the tri_mesh
		for( TriangleMesh::face_iter fit = tri_mesh->faces_begin(); fit != tri_mesh->faces_end(); ++fit )
		{
		  TriangularFace triface;
		  int i = 0;
		  for( TriangleMesh::FaceVertexIter fvit = tri_mesh->fv_iter(*fit); fvit; ++fvit, ++i )
		  {
		    assert( i >= 0 ); assert( i < 3 );
		    triface.idx[i] = getNumVerts()+(*fvit).idx();
		  }
		  m_faces.push_back(triface);
		  m_face_radii.push_back(0.0);
		}

		// Extract the vertex radii from the tri_mesh
		for( int i = 0; i < tri_mesh->nv(); ++i ) m_vertex_radii.push_back(0.0);

		// Extract masses from the tri_mesh (just infinity for this object)
		for( int i = 0; i < tri_mesh->nv(); ++i ) m_masses.push_back(std::numeric_limits<double>::infinity());

		// Update vector that tracks the ScriptedTriangleMesh DOF in the system
		m_base_triangle_indices.push_back(getNumDof());
		
		// Update total number of DOF in the system
		m_num_dof += 3*tri_mesh->nv();

	}

  m_xn.resize(getNumDof()); m_xn.setConstant(std::numeric_limits<double>::signaling_NaN());
  m_xnp1.resize(getNumDof()); m_xnp1.setConstant(std::numeric_limits<double>::signaling_NaN());
  m_vnphalf.resize(getNumDof()); m_vnphalf.setConstant(std::numeric_limits<double>::signaling_NaN());
  
  prepareForExecution();

}

RodTimeStepper* BridsonStepper::getInternalRodTimeStepper(int i) {
	assert( m_steppers.size() > i );
	
	return m_steppers[i];
}


void BridsonStepper::addTriangleMesh( TriangleMesh* tri_mesh )
{
  assert( tri_mesh != NULL );
  assert( (int) m_vertex_radii.size() == getNumVerts() );

  
  // Extract faces from the tri_mesh
  for( TriangleMesh::face_iter fit = tri_mesh->faces_begin(); fit != tri_mesh->faces_end(); ++fit )
  {
    TriangularFace triface;
    int i = 0;
    for( TriangleMesh::FaceVertexIter fvit = tri_mesh->fv_iter(*fit); fvit; ++fvit, ++i )
    {
      assert( i >= 0 ); assert( i < 3 );
      triface.idx[i] = getNumVerts()+(*fvit).idx();
    }
    m_faces.push_back(triface);
    m_face_radii.push_back(0.0);
  }
  assert( m_faces.size() == m_face_radii.size() );
  
  // Extract the vertex radii from the tri_mesh
  for( int i = 0; i < tri_mesh->nv(); ++i ) m_vertex_radii.push_back(0.0);

  // Extract masses from the tri_mesh (just infinity for this object)
  for( int i = 0; i < tri_mesh->nv(); ++i ) m_masses.push_back(std::numeric_limits<double>::infinity());

  // Add the mesh to the system
  m_triangle_meshes.push_back(tri_mesh);

  // Update vector that tracks the ScriptedTriangleMesh DOF in the system
  m_base_triangle_indices.push_back(getNumDof());
  
  // Update total number of DOF in the system
  m_num_dof += 3*tri_mesh->nv();

  // Resize the internal temporary storage
  m_xn.resize(getNumDof()); m_xn.setConstant(std::numeric_limits<double>::signaling_NaN());
  m_xnp1.resize(getNumDof()); m_xnp1.setConstant(std::numeric_limits<double>::signaling_NaN());
  m_vnphalf.resize(getNumDof()); m_vnphalf.setConstant(std::numeric_limits<double>::signaling_NaN());
  
  //m_bvh->regenerateBVH( m_xn, m_edges, m_faces, m_vertex_radii );
  //m_bvh->regenerateBVH();
  
  assert( m_base_triangle_indices.size() == m_triangle_meshes.size() );
  assert( (int) m_masses.size() == getNumVerts() );
  assert( (int) m_vertex_radii.size() == getNumVerts() );
}


void BridsonStepper::exertPenaltyImpulses( std::vector<EdgeEdgeProximityCollision>& edg_edg_cllsns, std::vector<VertexFaceProximityCollision>& vrtx_fce_cllsns, VecXd& v )
{
  // Exert edge-edge penalty (spring) impulses based on pre-step collisions
  for( int i = 0; i < (int) edg_edg_cllsns.size(); ++i )
  {
    assert( m_edg_edg_pnlty >= 0.0 );
    
    int idxa0 = edg_edg_cllsns[i].e0_v0;
    int idxa1 = edg_edg_cllsns[i].e0_v1;
    int idxb0 = edg_edg_cllsns[i].e1_v0;
    int idxb1 = edg_edg_cllsns[i].e1_v1;
    
    // Compute the relative velocity of the edges
    Vec3d relvel = computeRelativeVelocity( m_vnphalf, idxa0, idxa1, idxb0, idxb1, edg_edg_cllsns[i].s, edg_edg_cllsns[i].t );
    double relvelnorm = relvel.dot(edg_edg_cllsns[i].n);
    if( relvelnorm > 0.0 ) continue;
    
    double mass_sum_a = m_masses[idxa0]+m_masses[idxa1];
    // m_edg_edg_pnlty*edg_edg_cllsns[i].pen*m_dt
    double magnitude_a = std::min( m_edg_edg_pnlty*edg_edg_cllsns[i].pen*m_dt, std::numeric_limits<double>::infinity() + mass_sum_a*(0.25*edg_edg_cllsns[i].pen/m_dt - relvelnorm) );
    Vec3d impulse_a = -magnitude_a*edg_edg_cllsns[i].n;

    double mass_sum_b = m_masses[idxb0]+m_masses[idxb1];
    double magnitude_b = std::min( m_edg_edg_pnlty*edg_edg_cllsns[i].pen*m_dt, std::numeric_limits<double>::infinity() + mass_sum_b*(0.25*edg_edg_cllsns[i].pen/m_dt - relvelnorm) );
    Vec3d impulse_b = magnitude_b*edg_edg_cllsns[i].n;
    
    exertEdgeImpulse(  impulse_a, m_masses[idxa0], m_masses[idxa1], edg_edg_cllsns[i].s, idxa0, idxa1, m_vnphalf );
    exertEdgeImpulse(  impulse_b, m_masses[idxb0], m_masses[idxb1], edg_edg_cllsns[i].t, idxb0, idxb1, m_vnphalf );
  }
  
  // Exert vertex-face penalty (spring) impulses based on pre-step collisions
  for( int i = 0; i < (int) vrtx_fce_cllsns.size(); ++i )
  {
    //std::cout << "NumRods: " << getNumRods() << std::endl;
    //std::cout << "Vertex-face collisions: " << vrtx_fce_cllsns.size() << std::endl;
    //
    //for( int j = 0; j < (int) vrtx_fce_cllsns.size(); ++j )
    //{
    //  std::cout << "   cllsn: " << j << std::endl;
    //  std::cout << "       v0: " << vrtx_fce_cllsns[j].v0 << std::endl;
    //  std::cout << "       f0: " << vrtx_fce_cllsns[j].f0 << std::endl;
    //  std::cout << "       f1: " << vrtx_fce_cllsns[j].f1 << std::endl;
    //  std::cout << "       f2: " << vrtx_fce_cllsns[j].f2 << std::endl;
    //}
    
    assert( m_vrt_fc_pnlty >= 0.0 );
    double magnitude = m_vrt_fc_pnlty*vrtx_fce_cllsns[i].pen*m_dt;
    Vec3d impulse = magnitude*vrtx_fce_cllsns[i].n;
    
    int vrtidx = vrtx_fce_cllsns[i].v0;
    int fcidx0 = vrtx_fce_cllsns[i].f0;
    int fcidx1 = vrtx_fce_cllsns[i].f1;
    int fcidx2 = vrtx_fce_cllsns[i].f2;
    
    // Compute the relative velocity of the edges
    Vec3d relvel = computeRelativeVelocity( m_vnphalf, vrtidx, fcidx0, fcidx1, fcidx2, vrtx_fce_cllsns[i].u, vrtx_fce_cllsns[i].v, vrtx_fce_cllsns[i].w );
    
    //std::cout << "I: " << impulse.transpose() << std::endl;
    //std::cout << "relvel: " << relvel.transpose() << std::endl;
    
    if( relvel.dot(vrtx_fce_cllsns[i].n) <= 0.0 )
    {
      exertVertexImpulse( impulse, m_masses[vrtidx], vrtidx, m_vnphalf );      
	  	exertVertexFrictionImpulse( impulse, m_masses[vrtidx], vrtidx, m_vnphalf, relvel, vrtx_fce_cllsns[i].n );      

      exertFaceImpulse( -impulse, m_masses[fcidx0], m_masses[fcidx1], m_masses[fcidx2], 
                        vrtx_fce_cllsns[i].u, vrtx_fce_cllsns[i].v, vrtx_fce_cllsns[i].w, 
                        fcidx0, fcidx1, fcidx2, m_vnphalf );
    }
  }
}


void BridsonStepper::exertInelasticImpulses( std::vector<EdgeEdgeContinuousTimeCollision>& edg_edg_cllsns, std::vector<VertexFaceContinuousTimeCollision>& vrtx_fce_cllsns, VecXd& v )
{
  // Inelastic edge-edge impulses
  for( int i = 0; i < (int) edg_edg_cllsns.size(); ++i )
  {
    int idxa0 = edg_edg_cllsns[i].e0_v0;
    int idxa1 = edg_edg_cllsns[i].e0_v1;
    int idxb0 = edg_edg_cllsns[i].e1_v0;
    int idxb1 = edg_edg_cllsns[i].e1_v1;

    // Compute the relative velocity of the edges at the collision point
    Vec3d relvel = computeRelativeVelocity( m_vnphalf, idxa0, idxa1, idxb0, idxb1, edg_edg_cllsns[i].s, edg_edg_cllsns[i].t );
    double magrelvel = relvel.dot(edg_edg_cllsns[i].n);

    // If the edges are still approaching at the collision point
    if( magrelvel < 0.0 )
    {
      //std::cout << "n: " << edge_edge_collisions[i].n.transpose() << std::endl;
      Vec3d I = computeEdgeEdgeInelasticImpulse( m_masses[idxa0], m_masses[idxa1], m_masses[idxb0], m_masses[idxb1],
                                                 edg_edg_cllsns[i].s, edg_edg_cllsns[i].t, magrelvel, edg_edg_cllsns[i].n );

      exertEdgeImpulse( -I, m_masses[idxa0], m_masses[idxa1], edg_edg_cllsns[i].s, idxa0, idxa1, m_vnphalf );
      exertEdgeImpulse(  I, m_masses[idxb0], m_masses[idxb1], edg_edg_cllsns[i].t, idxb0, idxb1, m_vnphalf );
    }
  }

  // Vertex-face inelastic impulses
  for( int i = 0; i < (int) vrtx_fce_cllsns.size(); ++i )
  {
    int vrt_idx  = vrtx_fce_cllsns[i].v0;
    int fce_idx0 = vrtx_fce_cllsns[i].f0;
    int fce_idx1 = vrtx_fce_cllsns[i].f1;
    int fce_idx2 = vrtx_fce_cllsns[i].f2;

    double u = vrtx_fce_cllsns[i].u;
    double v = vrtx_fce_cllsns[i].v;
    double w = vrtx_fce_cllsns[i].w;
    
    // Compute the relative velocity of the edges at the collision point
    Vec3d relvel = computeRelativeVelocity(m_vnphalf,vrt_idx,fce_idx0,fce_idx1,fce_idx2,u,v,w);
    double magrelvel = relvel.dot(vrtx_fce_cllsns[i].n);

    if( magrelvel < 0.0 )
    {
      //Vec3d I = computeVertexFaceInelasticImpulse( m_masses[vrt_idx], m_masses[fce_idx0], m_masses[fce_idx1], m_masses[fce_idx2],
//                                                   u, v, w, magrelvel, vrtx_fce_cllsns[i].n );

	  Vec3d I = -magrelvel * vrtx_fce_cllsns[i].n * m_masses[vrt_idx];

//      exertFaceImpulse( -I, m_masses[fce_idx0], m_masses[fce_idx1], m_masses[fce_idx2], 
//                        u, v, w, fce_idx0, fce_idx1, fce_idx2, m_vnphalf );
//      //std::cout << m_vnphalf.transpose() << std::endl;
//		std::cout << I << "\n";
		
      exertVertexImpulse( I, m_masses[vrt_idx], vrt_idx, m_vnphalf );
//      //std::cout << m_vnphalf.transpose() << std::endl;

			exertVertexFrictionImpulse( I, m_masses[vrt_idx], vrt_idx, m_vnphalf, relvel, vrtx_fce_cllsns[i].n );      

    }
  }  
}

// For MIT's sewing machine project
// Only ONE rod
// It will solve dynamic of rods for the same time step "twice".
// First, run without boundary condition preventing a vertex from falling down under the belt.
// Then, it adds boundary condition for each vertex ended up with a position under the belt.
// Restore all configurations and run the same time step with boundary conditions.
void BridsonStepper::execute2()
{
	;
}


void BridsonStepper::execute()
{
//	std::cout << "\nNew Step\n";
	if (m_ground_fix && (int)m_triangle_meshes.size() > 0) 
	{
		execute2();
		return;
	}

	std::cout << "\nex 1\n";
	
  assert( m_edges.size() == m_edge_radii.size() );
  assert( (int) m_masses.size() == m_xn.size()/3 );
  assert( m_xn.size() == m_xnp1.size() );
  assert( m_xn.size() == m_vnphalf.size() );
  
  if( m_respns_enbld )
  {
    // Sanity check to ensure rods are not "internally colliding"
    for( int i = 0; i < (int) m_rods.size(); ++i ) ensureNoCollisionsByDefault( *m_rods[i] );

    // Pre time step position
    // TODO: by rearranging some copies, we can eliminate this extractPositions call.
    extractPositions(m_rods,m_base_indices,m_xn);
  }
    
//	std::cout << "\n\n\nstart\n" << m_rods[0]->getTheta(0) << " " << m_rods[0]->getTheta(1) << "\n";
    
  // Step rods forward ignoring collisions
  START_TIMER("BridsonStepperDynamics");

  for( int i = 0; i < (int) m_steppers.size(); ++i ) m_steppers[i]->execute();

  STOP_TIMER("BridsonStepperDynamics")

//	std::cout << "after dynamic : \n" << m_rods[0]->getTheta(0) << " " << m_rods[0]->getTheta(1) << "\n";

  // Move objects (mesh)
  for( int i = 0; i < (int) m_triangle_meshes.size(); ++i )
  {
    for( TriangleMesh::vertex_iter vit = m_triangle_meshes[i]->vertices_begin(); vit != m_triangle_meshes[i]->vertices_end(); ++vit )
    {
      m_triangle_meshes[i]->getVertex(*vit) += m_dt * m_triangle_meshes[i]->getVelocity(*vit);
    }
  }
  

 	if (m_fix_material_frame) {
		for( int i = 0; i < (int) m_rods.size(); ++i ) 
		{
			bool lastVertFix = false;
			for( int j = 0; j < m_rods[i]->nv(); ++j )
			{
				Vec3d x = m_rods[i]->getVertex(j);
				double r = m_rods[i]->getRadiusScale() * m_rods[i]->radiusA(std::min(m_rods[i]->ne() - 1, j));
				
				if (x[1] < r) {
					
					if (lastVertFix && 0 < j && j < m_rods[i]->nv() - 1) {
						if (!m_rods[i]->getBoundaryCondition()->isEdgeScripted(j-1)) {
							m_rods[i]->getBoundaryCondition()->setDesiredEdgeAngle(j-1, m_rods[i]->getTheta(j-1));
						}
					}
					
					lastVertFix = true;
				} else {
					lastVertFix = false;
				}
			}
		}
	}
	  
  
  if( m_respns_enbld )
  {
    START_TIMER("BridsonStepperResponse");
    // Post time step position
    extractPositions(m_rods,m_base_indices,m_xnp1);
    
    // Average velocity over the timestep just completed
    m_vnphalf = (m_xnp1-m_xn)/m_dt;
    
    // TODO: The setup here is a little space inefficient (pass around a bunch of lists and stuff), rework it later
    if( 0 && m_pnlty_enbld )
    {
      // Detect possible "proximity" collisions based on pre-step positions
      std::vector<EdgeEdgeProximityCollision> edge_edge_collisions;
      std::vector<VertexFaceProximityCollision> vertex_face_collisions;
      m_bvh->getProximityCollisions(edge_edge_collisions,vertex_face_collisions);
      
      //if( vertex_face_collisions.size() != 0 ) std::cout << "Vertex Face: " << vertex_face_collisions.size() << std::endl;
      if( edge_edge_collisions.size() != 0 ) std::cout << "Edge Edge: " << edge_edge_collisions.size() << std::endl;
      
      // Determine if any of the possible collisions happen. For real yo.
      std::vector<EdgeEdgeProximityCollision> edge_edge_for_real;
      detectEdgeEdgeProximityCollisions( m_xn, edge_edge_collisions, edge_edge_for_real );
      
      std::vector<VertexFaceProximityCollision> vertex_face_for_real;
      detectVertexFaceProximityCollisions( m_xn, vertex_face_collisions, vertex_face_for_real );
      
      // Apply penalty impulses
      exertPenaltyImpulses( edge_edge_for_real, vertex_face_for_real, m_vnphalf );
    }
    
    if( m_itrv_inlstc_enbld && m_num_inlstc_itrns != 0 )
    {
      // Detect possible continuous-time collisions
      std::vector<EdgeEdgeContinuousTimeCollision> edge_edge_collisions;
      std::vector<VertexFaceContinuousTimeCollision> vertex_face_collisions;
      m_bvh->getContinuousTimeCollisions(edge_edge_collisions,vertex_face_collisions);
      
      //detectEdgeEdgeContinuousTimeCollisions( m_xn, m_vnphalf, edge_edge_collisions );
      //detectVertexFaceContinuousTimeCollisions( m_xn, m_vnphalf, vertex_face_collisions );

      // Determine if any of the possible collisions happen. For real yo.
      std::vector<EdgeEdgeContinuousTimeCollision> edge_edge_for_real;
//      detectEdgeEdgeContinuousTimeCollisions( m_xn, m_vnphalf, edge_edge_collisions, edge_edge_for_real );
      
      std::vector<VertexFaceContinuousTimeCollision> vertex_face_for_real;
      detectVertexFaceContinuousTimeCollisions( m_xn, m_vnphalf, vertex_face_collisions, vertex_face_for_real );
      
      // Iterativly apply inelastic impulses
      int itr = 0;
      for( itr = 0; itr < m_num_inlstc_itrns; ++itr )
      {
        if( edge_edge_for_real.size()+vertex_face_for_real.size() == 0 ) break;
        
        // TODO: Add a little extra kick here
        exertInelasticImpulses( edge_edge_for_real, vertex_face_for_real, m_vnphalf );
        
        edge_edge_collisions.clear();
        vertex_face_collisions.clear();
        m_bvh->getContinuousTimeCollisions(edge_edge_collisions,vertex_face_collisions);
        
        edge_edge_for_real.clear();
//        detectEdgeEdgeContinuousTimeCollisions( m_xn, m_vnphalf, edge_edge_collisions, edge_edge_for_real );
        vertex_face_for_real.clear();
        detectVertexFaceContinuousTimeCollisions( m_xn, m_vnphalf, vertex_face_collisions, vertex_face_for_real );
      }
      if( itr == m_num_inlstc_itrns ) std::cout << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Exceeded maximum number of inelastic iterations " << m_num_inlstc_itrns << std::endl;
    }
    
    // Compute final positions from corrected velocities
    m_xnp1 = m_xn + m_dt*m_vnphalf;
    
    // Copy new positions and velocities back to rods
    restorePositions( m_rods, m_xnp1 );
    restoreVelocities( m_rods, m_vnphalf );

    // Update frames and such in the rod (Is this correct? Will this do some extra stuff?)
    for( int i = 0; i < (int) m_rods.size(); ++i ) m_rods[i]->updateProperties();
    
    STOP_TIMER("BridsonStepperResponse");
  } // if( m_respns_enbld )
}
  
void BridsonStepper::enableReseponse()
{
  m_respns_enbld = true;
}

void BridsonStepper::disableResponse()
{
  m_respns_enbld = false;
}

void BridsonStepper::enablePenaltyImpulses() 
{ 
  m_pnlty_enbld = true; 
}


void BridsonStepper::disablePenaltyImpulses() 
{ 
  m_pnlty_enbld = false; 
}


void BridsonStepper::enableIterativeInelasticImpulses() 
{ 
  m_itrv_inlstc_enbld = true; 
}


void BridsonStepper::disableIterativeInelasticImpulses() 
{ 
  m_itrv_inlstc_enbld = false; 
}


void BridsonStepper::setNumInelasticIterations( const int& num_itr )
{
  assert( num_itr >= 0 );
  m_num_inlstc_itrns = num_itr;
}


void BridsonStepper::setEdgeEdgePenalty( const double& k )
{
  assert( k >= 0.0 );
  m_edg_edg_pnlty = k;
}


void BridsonStepper::setVertexFacePenalty( const double& k )
{
  assert( k >= 0.0 );
  m_vrt_fc_pnlty = k;
}


// Ensures each rod edge has circular cross section. 
void BridsonStepper::ensureCircularCrossSection( const ElasticRod& rod ) const
{
  // Ensure circular cross section
  for( int i = 0; i < (int) rod.ne(); ++i )
  {
    if( rod.getRadiusScale()*rod.radiusA(i) != rod.getRadiusScale()*rod.radiusB(i) )
    {
      std::cerr << "Contact currently not supported for non-circular cross sections. Exiting." << std::endl;
      assert( false );
    }
  }
}


// Ensures each internal rod edge has length less than sum of neighbors' radii. 
void BridsonStepper::ensureNoCollisionsByDefault( const ElasticRod& rod ) const
{
  // Ensure "non-attached" edges are not colliding by default
  for( int i = 1; i < (int) rod.ne()-1; ++i )
  {
    double edgelen = rod.getEdge(i).norm();
    double radsum = rod.getRadiusScale()*rod.radiusA(i-1) + rod.getRadiusScale()*rod.radiusA(i+1);
    if( edgelen <= radsum )
    {
      std::cerr << "Detected edges that collide by default. Exiting." << std::endl;
      assert( false );
    }
  }
}


int BridsonStepper::getNumDof() const
{
  assert( m_num_dof >= 0 );
  return m_num_dof;
}

  
int BridsonStepper::getNumVerts() const
{
  assert( m_num_dof%3 == 0 );
  return m_num_dof/3;  
}
  
  
RodTimeStepper* BridsonStepper::createRodTimeStepper( ElasticRod& rod )
{
  RodTimeStepper* stepper = new RodTimeStepper(rod);
  
  assert( m_method == RodTimeStepper::IMPL_EULER );
  assert( m_method == RodTimeStepper::SYMPL_EULER || m_method == RodTimeStepper::IMPL_EULER );
  stepper->setDiffEqSolver(m_method);
  
  assert( m_dt > 0.0 );
  stepper->setTimeStep(m_dt);
  
  if( m_mass_damping != 0.0 ) stepper->addExternalForce( new RodMassDamping(m_mass_damping) );
  if( m_gravity.norm() > 0.0 ) stepper->addExternalForce( new RodGravity(m_gravity) );
  
  if (m_ground_fix && m_friction_enbld && m_friction_force) {
  	stepper->addExternalForce( m_friction_force );
  }
  
  assert( m_max_implicit_iterations > 0 );
  stepper->setMaxIterations(m_max_implicit_iterations);
  
  return stepper;
}


void BridsonStepper::extractTheta( const std::vector<ElasticRod*>& rods, VecXd& angles, VecXd& angleDots )
{
	int size = 0;
	
  for( int i = 0; i < (int) rods.size(); ++i ) 
  	size += rods[i]->ne();
  	
  angles.resize(size);
  angleDots.resize(size);
  
  int c = 0;
  
  for( int i = 0; i < (int) rods.size(); ++i ) {
		for( int j = 0; j < rods[i]->ne(); ++j )
		{
			angles[c] = rods[i]->getTheta(j);
			angleDots[c] = rods[i]->getThetaDot(j);
			c++;
		}
	}
}


void BridsonStepper::extractPositions( const std::vector<ElasticRod*>& rods, const std::vector<int>& base_indices, VecXd& positions )
{
  //std::cout << "!!!!!!!!!!!!!!! EXTRACTING POSITIONS !!!!!!!!!!!!!!!" << std::endl;
  assert( rods.size() == base_indices.size() );
  assert( getNumDof() == positions.size() );

  if( getNumDof() == 0 ) return;
  
  // DEBUG: Fill with NaNs
  positions.setConstant(std::numeric_limits<double>::signaling_NaN());
  
  for( int i = 0; i < (int) rods.size(); ++i ) for( int j = 0; j < rods[i]->nv(); ++j )
  {
    assert( base_indices[i]+3*j+2 < positions.size() );
    positions.segment<3>(base_indices[i]+3*j) = rods[i]->getVertex(j);
  }
  
  assert( m_triangle_meshes.size() == m_base_triangle_indices.size() );
  //std::cout << "!!!!!!!!!!!!!!! EXTRACTING TRIANGLES !!!!!!!!!!!!!!!" << std::endl;

  for( int i = 0; i < (int) m_triangle_meshes.size(); ++i )
  {
    int j = 0;
    for( TriangleMesh::vertex_iter vit = m_triangle_meshes[i]->vertices_begin(); vit != m_triangle_meshes[i]->vertices_end(); ++vit, ++j )
    {
      assert( m_base_triangle_indices[i]+3*j+2 < positions.size() );
      positions.segment<3>(m_base_triangle_indices[i]+3*j) = m_triangle_meshes[i]->getVertex(*vit);
    }
  }
  
  assert( (positions.cwise()==positions).all() );
}


void BridsonStepper::extractVelocities( const std::vector<ElasticRod*>& rods, const std::vector<int>& base_indices, VecXd& velocities )
{
  assert( rods.size() == base_indices.size() );
  assert( getNumDof() == velocities.size() );
  
  if( getNumDof() == 0 ) return;
  
  // DEBUG: Fill with NaNs
  velocities.setConstant(std::numeric_limits<double>::signaling_NaN());
  
  for( int i = 0; i < (int) rods.size(); ++i ) for( int j = 0; j < rods[i]->nv(); ++j )
  {
    assert( base_indices[i]+3*j+2 < velocities.size() );
    velocities.segment<3>(base_indices[i]+3*j) = rods[i]->getVelocity(j);
  }
  
  assert( m_triangle_meshes.size() == m_base_triangle_indices.size() );
  for( int i = 0; i < (int) m_triangle_meshes.size(); ++i )
  {
    int j = 0;
    for( TriangleMesh::vertex_iter vit = m_triangle_meshes[i]->vertices_begin(); vit != m_triangle_meshes[i]->vertices_end(); ++vit, ++j )
    {
      assert( m_base_triangle_indices[i]+3*j+2 < velocities.size() );
      velocities.segment<3>(m_base_triangle_indices[i]+3*j) = Vec3d::Zero();
    }
  }

  assert( (velocities.cwise()==velocities).all() );
}


void BridsonStepper::restorePositions( std::vector<ElasticRod*>& rods, const VecXd& positions )
{
  assert( rods.size() == m_base_indices.size() );
  
  for( int i = 0; i < (int) m_base_indices.size(); ++i )
  {
    for( int j = 0; j < rods[i]->nv(); ++j )
    {
      rods[i]->setVertex(j,positions.segment<3>(m_base_indices[i]+3*j));
    }
  }
}


void BridsonStepper::restoreVelocities( std::vector<ElasticRod*>& rods, const VecXd& velocities )
{
  assert( rods.size() == m_base_indices.size() );
  
  for( int i = 0; i < (int) m_base_indices.size(); ++i )
  {
    for( int j = 0; j < rods[i]->nv(); ++j )
    {
      rods[i]->setVelocity(j,velocities.segment<3>(m_base_indices[i]+3*j));
    }
  }
}

void BridsonStepper::restoreTheta( const std::vector<ElasticRod*>& rods, const VecXd& angles, const VecXd& angleDots  )
{
  int c = 0;
  
  for( int i = 0; i < (int) rods.size(); ++i ) {
		for( int j = 0; j < rods[i]->ne(); ++j )
		{
			rods[i]->setTheta(j, angles[c]);
			rods[i]->setThetaDot(j, angleDots[c]);
			c++;
		}
	}
}


bool BridsonStepper::edgesShareVertex( const std::pair<int,int>& edgei, const std::pair<int,int>& edgej ) const
{
  assert( edgei.first < getNumVerts() );
  assert( edgei.second < getNumVerts() );
  assert( edgej.first < getNumVerts() );
  assert( edgej.second < getNumVerts() );

  return edgei.first == edgej.first || edgei.second == edgej.second || edgei.first == edgej.second || edgei.second == edgej.first;
}

bool BridsonStepper::edgesSharevertex( const int& e0v0, const int& e0v1, const int& e1v0, const int& e1v1 ) const
{
  return e0v0 == e1v0 || e0v0 == e1v1 || e0v1 == e1v0 || e0v1 == e1v1;
}
  

bool BridsonStepper::vertexAndFaceShareVertex( const int& vertex, const int& face ) const
{
  assert( face < (int) m_faces.size() );
  assert( vertex < getNumVerts() );
  assert( m_faces[face].idx[0] < getNumVerts() );
  assert( m_faces[face].idx[1] < getNumVerts() );
  assert( m_faces[face].idx[2] < getNumVerts() );

  if( m_faces[face].idx[0] == vertex ) return true;
  if( m_faces[face].idx[1] == vertex ) return true;
  if( m_faces[face].idx[2] == vertex ) return true;
  
  return false;
}

bool BridsonStepper::vertexAndFaceShareVertex( const int& v, const int& f0, const int& f1, const int& f2 ) const
{
  return v == f0 || v == f1 || v == f2;
}
  

// Generates a list of all possible edge-edge proximity colisions
void BridsonStepper::generateAllEdgeEdgeProximityCollisionPairs( std::vector<EdgeEdgeProximityCollision>& edge_edge_collisions ) const
{
  assert( false );
  
  assert( m_edges.size() == m_edge_radii.size() );
  assert( m_faces.size() == m_face_radii.size() );
  
  // All free-edge free-edge collisions
  for( int i = 0; i < (int) m_edges.size(); ++i ) for( int j = i+1; j < (int) m_edges.size(); ++j )
  {
    if( !edgesShareVertex(m_edges[i],m_edges[j]) ) 
    {
      EdgeEdgeProximityCollision pssbl_col;
      pssbl_col.e0_v0 = m_edges[i].first;
      pssbl_col.e0_v1 = m_edges[i].second;
      pssbl_col.e1_v0 = m_edges[j].first;
      pssbl_col.e1_v1 = m_edges[j].second;
      pssbl_col.r0 = m_edge_radii[i];
      pssbl_col.r1 = m_edge_radii[j];
      edge_edge_collisions.push_back(pssbl_col);
    }
  }
  
  // All face-edge free-edge collisions
  for( int i = 0; i < (int) m_faces.size(); ++i ) for( int j = 0; j < (int) m_edges.size(); ++j )
  {
    // First edge of the face
    EdgeEdgeProximityCollision pssbl_col;
    pssbl_col.e0_v0 = m_faces[i].idx[0];
    pssbl_col.e0_v1 = m_faces[i].idx[1];
    pssbl_col.e1_v0 = m_edges[j].first;
    pssbl_col.e1_v1 = m_edges[j].second;
    pssbl_col.r0 = m_face_radii[i];
    pssbl_col.r1 = m_edge_radii[j];
    edge_edge_collisions.push_back(pssbl_col);
   
    // Second edge of the face
    pssbl_col.e0_v0 = m_faces[i].idx[1];
    pssbl_col.e0_v1 = m_faces[i].idx[2];
    //pssbl_col.e1_v0 = m_edges[j].first;
    //pssbl_col.e1_v1 = m_edges[j].second;
    pssbl_col.r0 = m_face_radii[i];
    //pssbl_col.r1 = m_edge_radii[j];
    edge_edge_collisions.push_back(pssbl_col);

    // Third edge of the face
    pssbl_col.e0_v0 = m_faces[i].idx[2];
    pssbl_col.e0_v1 = m_faces[i].idx[0];
    //pssbl_col.e1_v0 = m_edges[j].first;
    //pssbl_col.e1_v1 = m_edges[j].second;
    pssbl_col.r0 = m_face_radii[i];
    //pssbl_col.r1 = m_edge_radii[j];
    edge_edge_collisions.push_back(pssbl_col);
  }
}

void BridsonStepper::generateAllEdgeEdgeContinuousTimeCollisionPairs( std::vector<EdgeEdgeContinuousTimeCollision>& edge_edge_collisions ) const
{
  assert( false );
  
  // All free-edge free-edge collisions
  for( int i = 0; i < (int) m_edges.size(); ++i ) for( int j = i+1; j < (int) m_edges.size(); ++j )
  {
    if( !edgesShareVertex(m_edges[i],m_edges[j]) ) 
    {
      EdgeEdgeContinuousTimeCollision pssbl_col;
      pssbl_col.e0_v0 = m_edges[i].first;
      pssbl_col.e0_v1 = m_edges[i].second;
      pssbl_col.e1_v0 = m_edges[j].first;
      pssbl_col.e1_v1 = m_edges[j].second;
      edge_edge_collisions.push_back(pssbl_col);
    }
  }
  
  // All face-edge free-edge collisions
  for( int i = 0; i < (int) m_faces.size(); ++i ) for( int j = 0; j < (int) m_edges.size(); ++j )
  {
    // First edge of the face
    EdgeEdgeContinuousTimeCollision pssbl_col;
    pssbl_col.e0_v0 = m_faces[i].idx[0];
    pssbl_col.e0_v1 = m_faces[i].idx[1];
    pssbl_col.e1_v0 = m_edges[j].first;
    pssbl_col.e1_v1 = m_edges[j].second;
    edge_edge_collisions.push_back(pssbl_col);
    
    // Second edge of the face
    pssbl_col.e0_v0 = m_faces[i].idx[1];
    pssbl_col.e0_v1 = m_faces[i].idx[2];
    //pssbl_col.e1_v0 = m_edges[j].first;
    //pssbl_col.e1_v1 = m_edges[j].second;
    edge_edge_collisions.push_back(pssbl_col);
    
    // Third edge of the face
    pssbl_col.e0_v0 = m_faces[i].idx[2];
    pssbl_col.e0_v1 = m_faces[i].idx[0];
    //pssbl_col.e1_v0 = m_edges[j].first;
    //pssbl_col.e1_v1 = m_edges[j].second;
    edge_edge_collisions.push_back(pssbl_col);
  }
}
  

void BridsonStepper::generateAllVertexFaceProximityCollisionPairs( std::vector<VertexFaceProximityCollision>& vertex_face_collisions ) const
{
  assert( false );
  
  for( int i = 0; i < getNumVerts(); ++i ) for( int j = 0; j < (int) m_faces.size(); ++j )
  {
    if( !vertexAndFaceShareVertex(i,j) ) 
    {
      VertexFaceProximityCollision pssbl_col;
      pssbl_col.v0 = i;
      pssbl_col.f0 = m_faces[j].idx[0];
      pssbl_col.f1 = m_faces[j].idx[1];
      pssbl_col.f2 = m_faces[j].idx[2];
      pssbl_col.r0 = m_vertex_radii[i];
      pssbl_col.r1 = m_face_radii[j];
      vertex_face_collisions.push_back( pssbl_col );
    }
  }
}


void BridsonStepper::generateAllVertexFaceContinuousTimeCollisionPairs( std::vector<VertexFaceContinuousTimeCollision>& vertex_face_collisions ) const
{
  assert( false );
  
  for( int i = 0; i < getNumVerts(); ++i ) for( int j = 0; j < (int) m_faces.size(); ++j )
  {
    if( !vertexAndFaceShareVertex(i,j) ) 
    {
      VertexFaceContinuousTimeCollision pssbl_col;
      pssbl_col.v0 = i;
      pssbl_col.f0 = m_faces[j].idx[0];
      pssbl_col.f1 = m_faces[j].idx[1];
      pssbl_col.f2 = m_faces[j].idx[2];
      vertex_face_collisions.push_back( pssbl_col );
    }
  }
}


void BridsonStepper::detectEdgeEdgeProximityCollisions( const VecXd& x, std::vector<EdgeEdgeProximityCollision>& pssbl_cllsns, std::vector<EdgeEdgeProximityCollision>& cllsns ) const
{
  assert( (int) m_masses.size() == getNumVerts() );
  assert( x.size() == getNumDof() );
  
  //std::vector<std::pair<int,int> > possible_edge_edge_collisions;
  //std::vector<EdgeEdgeProximityCollision> pssbl_cllsns;
  //generateAllEdgeEdgeProximityCollisionPairs( pssbl_cllsns );
  
  //std::cout << "Possible edge edge collisions: " << possible_edge_edge_collisions.size() << std::endl;
  
  for( int i = 0; i < (int) pssbl_cllsns.size(); ++i )
  {
    assert( pssbl_cllsns[i].e0_v0 < getNumVerts() ); assert( pssbl_cllsns[i].e0_v1 < getNumVerts() );
    assert( pssbl_cllsns[i].e1_v0 < getNumVerts() ); assert( pssbl_cllsns[i].e1_v1 < getNumVerts() );
    assert( pssbl_cllsns[i].r0 >= 0.0 ); assert( pssbl_cllsns[i].r1 >= 0.0 );

    int idxp1 = pssbl_cllsns[i].e0_v0;
    int idxq1 = pssbl_cllsns[i].e0_v1;
    int idxp2 = pssbl_cllsns[i].e1_v0;
    int idxq2 = pssbl_cllsns[i].e1_v1;

    if( edgesSharevertex(idxp1,idxq1,idxp2,idxq2) ) continue;
    
    // If both rods have a fixed vertex, ignore this particular collision
    if( (m_masses[idxp1]==std::numeric_limits<double>::infinity()||m_masses[idxq1]==std::numeric_limits<double>::infinity()) &&
        (m_masses[idxp2]==std::numeric_limits<double>::infinity()||m_masses[idxq2]==std::numeric_limits<double>::infinity())  )
    {
      continue;
    }

    const Vec3d p1 = x.segment<3>(3*idxp1);
    const Vec3d q1 = x.segment<3>(3*idxq1);
    const Vec3d p2 = x.segment<3>(3*idxp2);
    const Vec3d q2 = x.segment<3>(3*idxq2);
    
    Vec3d c1;
    Vec3d c2;    
    double sqrdist = ClosestPtSegmentSegment(p1,q1,p2,q2,pssbl_cllsns[i].s,pssbl_cllsns[i].t,c1,c2);
    
    assert( pssbl_cllsns[i].s == pssbl_cllsns[i].s ); assert( pssbl_cllsns[i].t == pssbl_cllsns[i].t );
    assert( pssbl_cllsns[i].s >= 0.0 ); assert( pssbl_cllsns[i].s <= 1.0 );
    assert( pssbl_cllsns[i].t >= 0.0 ); assert( pssbl_cllsns[i].t <= 1.0 );
    
    if( sqrdist < (pssbl_cllsns[i].r0+pssbl_cllsns[i].r1)*(pssbl_cllsns[i].r0+pssbl_cllsns[i].r1) )
    {
      pssbl_cllsns[i].pen = pssbl_cllsns[i].r0+pssbl_cllsns[i].r1-sqrt(sqrdist);      
      assert( pssbl_cllsns[i].pen > 0.0 );
      
      pssbl_cllsns[i].n = c2-c1;
      assert( pssbl_cllsns[i].n.norm() > 0.0 );
      
      pssbl_cllsns[i].n.normalize();
      assert( fabs(pssbl_cllsns[i].n.norm()-1.0) < 1.0e-6 );

      if( pssbl_cllsns[i].s != 0 && pssbl_cllsns[i].s != 1 && sqrt(sqrdist) > 1.0e-9 ) assert( fabs(pssbl_cllsns[i].n.dot(q1-p1)) < 1.0e-6 );
      if( pssbl_cllsns[i].t != 0 && pssbl_cllsns[i].t != 1 && sqrt(sqrdist) > 1.0e-9 ) assert( fabs(pssbl_cllsns[i].n.dot(q2-p2)) < 1.0e-6 );

      cllsns.push_back(pssbl_cllsns[i]);
    }
  }
}
  
  
void BridsonStepper::detectVertexFaceProximityCollisions( const VecXd& x, std::vector<VertexFaceProximityCollision>& pssbl_cllsns, std::vector<VertexFaceProximityCollision>& vetex_face_collisions ) const
{
  assert( (int) m_masses.size() == getNumVerts() );
  assert( x.size() == getNumDof() );
  assert( m_faces.size() == m_face_radii.size() );

  //std::vector<VertexFaceProximityCollision> pssbl_cllsns;
  //generateAllVertexFaceProximityCollisionPairs( pssbl_cllsns );
  
  for( int i = 0; i < (int) pssbl_cllsns.size(); ++i )
  {
    assert( pssbl_cllsns[i].v0 < getNumVerts() ); assert( pssbl_cllsns[i].f0 < getNumVerts() );
    assert( pssbl_cllsns[i].f1 < getNumVerts() ); assert( pssbl_cllsns[i].f2 < getNumVerts() );
    assert( pssbl_cllsns[i].r0 >= 0.0 ); assert( pssbl_cllsns[i].r1 >= 0.0 );

    int vrtidx = pssbl_cllsns[i].v0;
    int fcidx0 = pssbl_cllsns[i].f0;
    int fcidx1 = pssbl_cllsns[i].f1;
    int fcidx2 = pssbl_cllsns[i].f2;
    
    if( vertexAndFaceShareVertex( vrtidx, fcidx0, fcidx1, fcidx2 ) ) continue;

    // TODO: Add check for both having fixed vertices

    const Vec3d p1 = x.segment<3>(3*vrtidx);
    const Vec3d t0 = x.segment<3>(3*fcidx0);
    const Vec3d t1 = x.segment<3>(3*fcidx1);
    const Vec3d t2 = x.segment<3>(3*fcidx2);
    
    Vec3d cp = ClosestPtPointTriangle(p1,t0,t1,t2);
    double sqrdist = (p1-cp).norm();

    if( sqrdist < (pssbl_cllsns[i].r0+pssbl_cllsns[i].r1)*(pssbl_cllsns[i].r0+pssbl_cllsns[i].r1) )
    {
      Barycentric( t0, t1, t2, p1, pssbl_cllsns[i].u, pssbl_cllsns[i].v, pssbl_cllsns[i].w );

      assert( pssbl_cllsns[i].u >= 0.0 ); assert( pssbl_cllsns[i].u <= 1.0 );
      assert( pssbl_cllsns[i].v >= 0.0 ); assert( pssbl_cllsns[i].v <= 1.0 );
      assert( pssbl_cllsns[i].w >= 0.0 ); assert( pssbl_cllsns[i].w <= 1.0 );

      pssbl_cllsns[i].pen = pssbl_cllsns[i].r0+pssbl_cllsns[i].r1-sqrt(sqrdist);      
      assert( pssbl_cllsns[i].pen > 0.0 );

      pssbl_cllsns[i].n = p1-cp;
      assert( pssbl_cllsns[i].n.norm() > 0.0 );

      pssbl_cllsns[i].n.normalize();
      assert( fabs(pssbl_cllsns[i].n.norm()-1.0) < 1.0e-6 );

      // TODO: Add some checks that if the closest point is inside the triangle, normal is normal to each edge

      vetex_face_collisions.push_back(pssbl_cllsns[i]);      
    }
  }
}

  
// Detects all edge-edge continuous time collisions
void BridsonStepper::detectEdgeEdgeContinuousTimeCollisions( const VecXd& x, const VecXd& v, std::vector<EdgeEdgeContinuousTimeCollision>& pssbl_cllsns, std::vector<EdgeEdgeContinuousTimeCollision>& edge_edge_collisions )
{
  assert( x.size() == getNumDof() ); assert( v.size() == getNumDof() );
  assert( edge_edge_collisions.size() == 0 );
  
  //std::vector<EdgeEdgeContinuousTimeCollision> pssbl_cllsns;
  //generateAllEdgeEdgeContinuousTimeCollisionPairs( pssbl_cllsns );

  for( int i = 0; i < (int) pssbl_cllsns.size(); ++i )
  {
    assert( pssbl_cllsns[i].e0_v0 < getNumVerts() ); assert( pssbl_cllsns[i].e0_v1 < getNumVerts() );
    assert( pssbl_cllsns[i].e1_v0 < getNumVerts() ); assert( pssbl_cllsns[i].e1_v1 < getNumVerts() );

    int idxp1 = pssbl_cllsns[i].e0_v0;
    int idxq1 = pssbl_cllsns[i].e0_v1;
    int idxp2 = pssbl_cllsns[i].e1_v0;
    int idxq2 = pssbl_cllsns[i].e1_v1;
    
    if( edgesSharevertex(idxp1,idxq1,idxp2,idxq2) ) continue;
    
    // If both rods have a fixed vertex, ignore this particular collision
    if( (m_masses[idxp1]==std::numeric_limits<double>::infinity()||m_masses[idxq1]==std::numeric_limits<double>::infinity()) &&
        (m_masses[idxp2]==std::numeric_limits<double>::infinity()||m_masses[idxq2]==std::numeric_limits<double>::infinity())  )
    {
      continue;
    }

    const Vec3d p1 = x.segment<3>(3*idxp1);
    const Vec3d q1 = x.segment<3>(3*idxq1);
    const Vec3d p2 = x.segment<3>(3*idxp2);
    const Vec3d q2 = x.segment<3>(3*idxq2);
    
    const Vec3d vp1 = v.segment<3>(3*idxp1);
    const Vec3d vq1 = v.segment<3>(3*idxq1);
    const Vec3d vp2 = v.segment<3>(3*idxp2);
    const Vec3d vq2 = v.segment<3>(3*idxq2);

    std::vector<double> times;
    std::vector<double> errors;
    getCoplanarityTimes( p1, q1, p2, q2, p1+m_dt*vp1, q1+m_dt*vq1, p2+m_dt*vp2, q2+m_dt*vq2, times, errors );
    assert( times.size() == errors.size() );
    
    for( int j = 0; j < (int) times.size(); ++j )
    {
      if( times[j] != times[j] ) 
      { 
        if( !m_nan_enc )std::cout << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered NaN collision time from root finder. Supressing further messages of this type." << std::endl; 
        m_nan_enc = true;
        continue; 
      }
      if( times[j] == std::numeric_limits<double>::infinity() ) 
      { 
        if( !m_inf_enc ) std::cout << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered INF collision time from root finder. Supressing further messages of this type." << std::endl; 
        m_inf_enc = true;
        continue; 
      }
      if( times[j] < 0.0 ) 
      {
        if( !m_lt0_enc ) std::cout << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered scaled collision time " << times[j] << " less than 0.0. Supressing further messages of this type." << std::endl; 
        m_lt0_enc = true;
        continue; 
      }
      if( times[j] > 1.0 ) 
      {
        if( !m_gt0_enc ) std::cout << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered scaled collision time " << times[j] << " greater than 1.0. Supressing further messages of this type." << std::endl; 
        m_gt0_enc = true;
        continue; 
      }
      //assert( times[j] >= 0.0 ); assert( times[j] <= 1.0 );
      
      double dtime = times[j]*m_dt;
      
      // Determine if the collision actually happens
      const Vec3d p1col = p1+dtime*vp1;
      const Vec3d q1col = q1+dtime*vq1;
      const Vec3d p2col = p2+dtime*vp2;
      const Vec3d q2col = q2+dtime*vq2;
      
      Vec3d c1;
      Vec3d c2;
      double sqrdist = ClosestPtSegmentSegment( p1col, q1col, p2col, q2col, pssbl_cllsns[i].s, pssbl_cllsns[i].t, c1, c2 );
      assert( sqrdist >= 0.0 );
      assert( pssbl_cllsns[i].s >= 0.0 ); assert( pssbl_cllsns[i].s <= 1.0 );
      assert( pssbl_cllsns[i].t >= 0.0 ); assert( pssbl_cllsns[i].t <= 1.0 );

      // If, when they are coplanar, the objects are sufficiently close, register a collision
      if( sqrt(sqrdist) < 1.0e-6 )
      {
        // Compute a collision normal at the time of the collision. For a first attempt,
        // take the cross product of the edges.
        pssbl_cllsns[i].n = (q1col-p1col).cross(q2col-p2col);
        
        // If the edges happen to be parallel
        if( pssbl_cllsns[i].n.norm() == 0.0 )
        {
          // Use the pre-timestep positions of the collision points to generate a collision normal
          pssbl_cllsns[i].n = ((1.0-pssbl_cllsns[i].t)*p2+pssbl_cllsns[i].t*q2)-((1.0-pssbl_cllsns[i].s)*p1+pssbl_cllsns[i].s*q1);
          
          // Crazy corner case! Handle it if it comes up :)
          if( pssbl_cllsns[i].n.norm() == 0.0 )
          {
            bool lazy_person_implemented_this_normal_computation = false;
            if( !lazy_person_implemented_this_normal_computation ) std::cout << "YOU HIT AN UNSOPPORTED CODE PATH. REALLY DEGENERATE." << std::endl;
            assert( lazy_person_implemented_this_normal_computation );
          }
        }
        
        pssbl_cllsns[i].n.normalize();
        assert( fabs(pssbl_cllsns[i].n.norm()-1.0) < 1.0e-6 );

        // Make sure the normal points from first edge to second edge
        if( pssbl_cllsns[i].n.dot( ((1.0-pssbl_cllsns[i].t)*p2+pssbl_cllsns[i].t*q2)-((1.0-pssbl_cllsns[i].s)*p1+pssbl_cllsns[i].s*q1) ) < 0.0 ) pssbl_cllsns[i].n *= -1.0;

        edge_edge_collisions.push_back(pssbl_cllsns[i]);
      }
    }
  }
  
  assert( edge_edge_collisions.size() <= pssbl_cllsns.size() );
}

void BridsonStepper::detectVertexFaceContinuousTimeCollisions( const VecXd& x, const VecXd& v, std::vector<VertexFaceContinuousTimeCollision>& pssbl_cllsns, std::vector<VertexFaceContinuousTimeCollision>& vertex_face_collisions )
{
  assert( x.size() == getNumDof() ); assert( v.size() == getNumDof() );
  assert( vertex_face_collisions.size() == 0 );
  
  //std::vector<VertexFaceContinuousTimeCollision> pssbl_cllsns;
  //generateAllVertexFaceContinuousTimeCollisionPairs( pssbl_cllsns );
  
  for( int i = 0; i < (int) pssbl_cllsns.size(); ++i )
  {
    assert( pssbl_cllsns[i].v0 < getNumVerts() ); assert( pssbl_cllsns[i].f0 < getNumVerts() );
    assert( pssbl_cllsns[i].f1 < getNumVerts() ); assert( pssbl_cllsns[i].f2 < getNumVerts() );
    assert( pssbl_cllsns[i].v0 != pssbl_cllsns[i].f0 ); assert( pssbl_cllsns[i].v0 != pssbl_cllsns[i].f1 ); assert( pssbl_cllsns[i].v0 != pssbl_cllsns[i].f2 );

    int vrtidx = pssbl_cllsns[i].v0;
    int fcidx0 = pssbl_cllsns[i].f0;
    int fcidx1 = pssbl_cllsns[i].f1;
    int fcidx2 = pssbl_cllsns[i].f2;

    if( vertexAndFaceShareVertex( vrtidx, fcidx0, fcidx1, fcidx2 ) ) continue;
    
    assert( vrtidx != fcidx0 ); assert( vrtidx != fcidx1 ); assert( vrtidx != fcidx2 );
    
    // TODO: Add checks for both having fixed vertices

    const Vec3d p  = x.segment<3>(3*vrtidx);
    const Vec3d f0 = x.segment<3>(3*fcidx0);
    const Vec3d f1 = x.segment<3>(3*fcidx1);
    const Vec3d f2 = x.segment<3>(3*fcidx2);
    
    const Vec3d vp  = v.segment<3>(3*vrtidx);
    const Vec3d vf0 = v.segment<3>(3*fcidx0);
    const Vec3d vf1 = v.segment<3>(3*fcidx1);
    const Vec3d vf2 = v.segment<3>(3*fcidx2);

    std::vector<double> times;
    std::vector<double> errors;
    getCoplanarityTimes( p, f0, f1, f2, p+m_dt*vp, f0+m_dt*vf0, f1+m_dt*vf1, f2+m_dt*vf2, times, errors );
    assert( times.size() == errors.size() );

    for( int j = 0; j < (int) times.size(); ++j )
    {
      if( times[j] != times[j] ) 
      { 
        if( !m_nan_enc )std::cout << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered NaN collision time from root finder. Supressing further messages of this type." << std::endl; 
        m_nan_enc = true;
        continue; 
      }
      if( times[j] == std::numeric_limits<double>::infinity() ) 
      { 
        if( !m_inf_enc ) std::cout << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered INF collision time from root finder. Supressing further messages of this type." << std::endl; 
        m_inf_enc = true;
        continue; 
      }
      if( times[j] < 0.0 ) 
      {
        if( !m_lt0_enc ) std::cout << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered scaled collision time " << times[j] << " less than 0.0. Supressing further messages of this type." << std::endl; 
        m_lt0_enc = true;
        continue; 
      }
      if( times[j] > 1.0 ) 
      {
        if( !m_gt0_enc ) std::cout << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered scaled collision time " << times[j] << " greater than 1.0. Supressing further messages of this type." << std::endl; 
        m_gt0_enc = true;
        continue; 
      }

      double dtime = times[j]*m_dt;
      
      // Determine if the collision actually happens
      const Vec3d pcol  = p +dtime*vp;
      const Vec3d f0col = f0+dtime*vf0;
      const Vec3d f1col = f1+dtime*vf1;
      const Vec3d f2col = f2+dtime*vf2;
            
      Vec3d cp = ClosestPtPointTriangle(pcol,f0col,f1col,f2col);
      
      double dist = sqrt((pcol-cp).norm());

      // If, when they are coplanar, the objects are sufficiently close, register a collision
      if( dist < 1.0e-6 )
      {        
        // Compute a collision normal at the time of the collision. For a first attempt,
        // take the cross product of two face edges.
        pssbl_cllsns[i].n = (f2col-f0col).cross(f1col-f0col);

        // If the normal is of magnitude 0 (triangle degenerates into a line, I believe)
        if( pssbl_cllsns[i].n.norm() == 0.0 )
        {
          bool lazy_person_implemented_this_normal_computation = false;
          if( !lazy_person_implemented_this_normal_computation ) std::cout << "YOU HIT AN UNSOPPORTED CODE PATH. REALLY DEGENERATE." << std::endl;
          assert( lazy_person_implemented_this_normal_computation );
        }

        pssbl_cllsns[i].n.normalize();
        assert( fabs(pssbl_cllsns[i].n.norm()-1.0) < 1.0e-6 );

        // Make sure the normal points from the face to the point
        if( pssbl_cllsns[i].n.dot( p-f0 ) < 0.0 ) pssbl_cllsns[i].n *= -1.0;
        
        Barycentric( f0col, f1col, f2col, pcol, pssbl_cllsns[i].u, pssbl_cllsns[i].v, pssbl_cllsns[i].w );
        assert( pssbl_cllsns[i].u >= 0.0 ); assert( pssbl_cllsns[i].u <= 1.0 );
        assert( pssbl_cllsns[i].v >= 0.0 ); assert( pssbl_cllsns[i].v <= 1.0 );
        assert( pssbl_cllsns[i].w >= 0.0 ); assert( pssbl_cllsns[i].w <= 1.0 );

        vertex_face_collisions.push_back(pssbl_cllsns[i]);
      }      
    }
  }
  
  assert( vertex_face_collisions.size() <= pssbl_cllsns.size() );
}
  
Vec3d BridsonStepper::computeRelativeVelocity( const VecXd& v, const int& idxa0, const int& idxa1, const int& idxb0, const int& idxb1, const double& s, const double& t )
{
  assert( v.size() == getNumDof() );
  assert( idxa0 < getNumVerts() ); assert( idxa1 < getNumVerts() );
  assert( idxb0 < getNumVerts() ); assert( idxb1 < getNumVerts() );
  assert( s >= 0.0 ); assert( s <= 1.0 );
  assert( t >= 0.0 ); assert( t <= 1.0 );
  
  const Vec3d v0 = v.segment<3>(3*idxa0);
  const Vec3d v1 = v.segment<3>(3*idxa1);
  const Vec3d v2 = v.segment<3>(3*idxb0);
  const Vec3d v3 = v.segment<3>(3*idxb1);
  
  return ((1.0-t)*v2+t*v3)-((1.0-s)*v0+s*v1);
}

  
Vec3d BridsonStepper::computeRelativeVelocity( const VecXd& vel, const int& vrtidx, const int& fcidx0, const int& fcidx1, const int& fcidx2, const double& u, const double& v, const double& w )
{
  assert( vel.size() == getNumDof() );
  assert( vrtidx < getNumVerts() );
  assert( fcidx0 < getNumVerts() );
  assert( fcidx1 < getNumVerts() );
  assert( fcidx2 < getNumVerts() );
  assert( u >= 0.0 ); assert( u <= 1.0 );
  assert( v >= 0.0 ); assert( v <= 1.0 );
  assert( w >= 0.0 ); assert( w <= 1.0 );
  
  const Vec3d vp  = vel.segment<3>(3*vrtidx);
  const Vec3d vt0 = vel.segment<3>(3*fcidx0);
  const Vec3d vt1 = vel.segment<3>(3*fcidx1);
  const Vec3d vt2 = vel.segment<3>(3*fcidx2);
  
  return vp - (u*vt0+v*vt1+w*vt2);
}
  
  
Vec3d BridsonStepper::computeEdgeEdgeInelasticImpulse( const double& ma0, const double& ma1, const double& mb0, const double& mb1,
                                                       const double& s, const double& t, const double& relvel, const Vec3d& n )
{
  assert( ma0 > 0.0 ); assert( ma1 > 0.0 ); assert( mb0 > 0.0 ); assert( mb1 > 0.0 );
  assert( s >= 0.0 ); assert( s <= 1.0 ); assert( t >= 0.0 ); assert( t <= 1.0 );
  assert( relvel <= 0.0 );
  
  // Assumes negative relative velocity
  Vec3d numerator = -relvel*n;
  double denominator = (1-s)*(1-s)/ma0 + s*s/ma1 + (1-t)*(1-t)/mb0 + t*t/mb1;
  assert( denominator != 0.0 );
  return numerator/denominator;
}

  
Vec3d BridsonStepper::computeVertexFaceInelasticImpulse( const double& mvrt, const double& mfc0, const double& mfc1, const double& mfc2,
                                                         const double& u, const double& v, const double& w, const double& relvel, const Vec3d& n )
{
  assert( mvrt >= 0.0 ); assert( mfc0 >= 0.0 ); assert( mfc1 >= 0.0 ); assert( mfc2 >= 0.0 );
  assert( u >= 0.0 ); assert( u <= 1.0 ); assert( v >= 0.0 ); assert( v <= 1.0 ); assert( w >= 0.0 ); assert( w <= 1.0 );
  assert( relvel <= 0.0 );
  
  Vec3d numerator = -relvel*n;
  double denominator = 1/mvrt + u*u/mfc0 + v*v/mfc1 + w*w/mfc2;
  assert( denominator != 0.0 );
  return numerator/denominator;
}
			
void BridsonStepper::exertVertexFrictionImpulse( const Vec3d& I, const double& m, const int& idx, VecXd& v, Vec3d& vRel, Vec3d& normal) {
	const double epsilon = 1.0e-6;  
  
  if (m_friction_enbld && I.norm() > epsilon) {
//  	Vec3d vRel = v.segment<3>(3*idx) - v_face;
    Vec3d vRelNormal = vRel.dot(normal) * normal;
    Vec3d vRelTangent = vRel - vRelNormal;  	

    Scalar tangentVel = vRelTangent.norm();
		if (tangentVel < epsilon) return;

		Scalar frictionMag = m_friction_cof * I.norm() / m;
		
		if (frictionMag > tangentVel) frictionMag = tangentVel;		// Make sure that a vertex will not move in the opposite direction.

		v.segment<3>(3*idx) += -frictionMag * (vRelTangent / tangentVel);
  } 

}
			  
void BridsonStepper::exertVertexImpulse( const Vec3d& I, const double& m, const int& idx, VecXd& v )
{
  assert( m > 0.0 );
  assert( idx < getNumVerts() );
  assert( v.size() == getNumDof() );
  
  v.segment<3>(3*idx) += I/m;
}
  
  
void BridsonStepper::exertEdgeImpulse( const Vec3d& I, const double& m0, const double& m1, const double& alpha, const int& idx0, const int& idx1, VecXd& v )
{
  assert( m0 > 0.0 ); assert( m1 > 0.0 );
  assert( alpha >= 0.0 ); assert( alpha <= 1.0 );
  assert( idx0 < getNumVerts() ); assert( idx1 < getNumVerts() );
  assert( v.size() == getNumDof() );
  
  v.segment<3>(3*idx0) += ((1-alpha)/m0)*I;
  v.segment<3>(3*idx1) += (alpha/m1)*I;
}

  
void BridsonStepper::exertFaceImpulse( const Vec3d& I, const double& m0, const double& m1, const double& m2, 
                                       const double& u, const double& v, const double& w, 
                                       const int& idx0, const int& idx1, const int& idx2, VecXd& vel )
{
  assert( m0 > 0.0 ); assert( m1 > 0.0 ); assert( m2 > 0.0 );
  assert( u >= 0.0 ); assert( u <= 1.0 );
  assert( v >= 0.0 ); assert( v <= 1.0 );
  assert( w >= 0.0 ); assert( w <= 1.0 );
  assert( idx0 < getNumVerts() ); assert( idx1 < getNumVerts() ); assert( idx2 < getNumVerts() );
  assert( vel.size() == getNumDof() );

  vel.segment<3>(3*idx0) += u*I/m0;
  vel.segment<3>(3*idx1) += v*I/m1;
  vel.segment<3>(3*idx2) += w*I/m2;
}


//bool CandidateCollision::getContinuousTimeEdgeEdge( const Vec3d &x00, const Vec3d &x01, const Vec3d &x10, const Vec3d &x11,
//                                                    const Vec3d &v00, const Vec3d &v01, const Vec3d &v10, const Vec3d &v11,
//                                                    double dt, Collision &collision ) const
//{
//  double epsilon = 1.0e-6;
//  
//  std::vector<Real> times, errors;
//  getCoplanarityTimes(x00, x01, x10, x11, x00+v00*dt, x01+v01*dt, x10+v10*dt, x11+v11*dt, times, errors);
//  
//  for( int a = 0; a < (int) times.size(); ++a )
//  {
//    // Coplanarity test code returns times scaled to range [0,1]
//    double t = times[a] * dt;
//    
//    Vec3d xt00 = x00 + t * v00;
//    Vec3d xt01 = x01 + t * v01;
//    Vec3d xt10 = x10 + t * v10;
//    Vec3d xt11 = x11 + t * v11;
//    
//    double s1, s2;
//    Vec3d c1, c2;
//    double distsqr = ClosestPtSegmentSegment( xt00, xt01, xt10, xt11, s1, s2, c1, c2 );
//
//    if( distsqr < epsilon*epsilon ) continue;
//
//    Vec3d normal = c2 - c1;
//    assert( normal.norm() > 0.0 );
//    normal.normalize();
//
//    return true;
//  }
//  
//  return false;
//}
  
  

// Adapted from some code on Robert Bridson's website, I believe

void BridsonStepper::addUnique( std::vector<double>& a, double e ) const
{
  for(unsigned int i=0; i<a.size(); ++i) if(a[i]==e) return;
  a.push_back(e);
}
  
double BridsonStepper::triple( const Vec3d& a, const Vec3d& b, const Vec3d& c ) const
{ 
  return a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0]); 
}
  
double BridsonStepper::signed_volume( const Vec3d& x0, const Vec3d& x1, const Vec3d& x2, const Vec3d& x3 ) const
{
  // Equivalent to triple(x1-x0, x2-x0, x3-x0), six times the signed volume of the tetrahedron.
  // But, for robustness, we want the result (up to sign) to be independent of the ordering.
  // And want it as accurate as possible...
  // But all that stuff is hard, so let's just use the common assumption that all coordinates are >0,
  // and do something reasonably accurate in fp.
  
  // This formula does almost four times too much multiplication, but if the coordinates are non-negative
  // it suffers in a minimal way from cancellation error.
  return ( x0[0]*(x1[1]*x3[2]+x3[1]*x2[2]+x2[1]*x1[2])
          +x1[0]*(x2[1]*x3[2]+x3[1]*x0[2]+x0[1]*x2[2])
          +x2[0]*(x3[1]*x1[2]+x1[1]*x0[2]+x0[1]*x3[2])
          +x3[0]*(x1[1]*x2[2]+x2[1]*x0[2]+x0[1]*x1[2]) )
  
  - ( x0[0]*(x2[1]*x3[2]+x3[1]*x1[2]+x1[1]*x2[2])
     +x1[0]*(x3[1]*x2[2]+x2[1]*x0[2]+x0[1]*x3[2])
     +x2[0]*(x1[1]*x3[2]+x3[1]*x0[2]+x0[1]*x1[2])
     +x3[0]*(x2[1]*x1[2]+x1[1]*x0[2]+x0[1]*x2[2]) );
}

// All roots returned in interval [0,1]. Assumed geometry followed a linear
// trajectory between x and xnew. 
void BridsonStepper::getCoplanarityTimes( const Vec3d& x0, const Vec3d& x1, const Vec3d& x2, const Vec3d& x3,
                                          const Vec3d& xnew0, const Vec3d& xnew1, const Vec3d& xnew2, const Vec3d& xnew3,
                                          std::vector<double>& times, std::vector<double>& errors ) const
{
  const double tol = 1e-8;
  times.clear();
  errors.clear();
  
  // cubic coefficients, A*t^3+B*t^2+C*t+D (for t in [0,1])
  Vec3d x03=x0-x3;
  Vec3d x13=x1-x3;
  Vec3d x23=x2-x3;
  Vec3d v03=(xnew0-xnew3)-x03;
  Vec3d v13=(xnew1-xnew3)-x13;
  Vec3d v23=(xnew2-xnew3)-x23;
  
  double A = triple(v03,v13,v23);
  double B = triple(x03,v13,v23)+triple(v03,x13,v23)+triple(v03,v13,x23);
  double C = triple(x03,x13,v23)+triple(x03,v13,x23)+triple(v03,x13,x23);
  double D = triple(x03,x13,x23);
  
  const double convergence_tol = tol*(std::fabs(A)+std::fabs(B)+std::fabs(C)+std::fabs(D));
  
  // find intervals to check, or just solve it if it reduces to a quadratic =============================
  std::vector<double> interval_times;
  double discriminant=B*B-3*A*C; // of derivative of cubic, 3*A*t^2+2*B*t+C, divided by 4 for convenience
  if(discriminant<=0){ // monotone cubic: only one root in [0,1] possible
    // so we just 
    interval_times.push_back(0);
    interval_times.push_back(1);
  }
  else
  { // positive discriminant, B!=0
    if(A==0)
    { // the cubic is just a quadratic, B*t^2+C*t+D ========================================
      discriminant=C*C-4*B*D; // of the quadratic
      if(discriminant<=0)
      {
        double t=-C/(2*B);
        if(t>=-tol && t<=1+tol)
        {
          t=clamp(t, 0., 1.);
          double val = std::fabs(signed_volume((1-t)*x0+t*xnew0,
                                             (1-t)*x1+t*xnew1,
                                             (1-t)*x2+t*xnew2,
                                             (1-t)*x3+t*xnew3));
          if (val < convergence_tol)
          {
            times.push_back(t);
          }
        }
      }
      else
      { // two separate real roots
        double t0, t1;
        if(C>0) t0=(-C-std::sqrt(discriminant))/(2*B);
        else    t0=(-C+std::sqrt(discriminant))/(2*B);
        t1=D/(B*t0);
        if(t1<t0) std::swap(t0,t1);
        if(t0>=-tol && t0<=1+tol)
        {
          times.push_back(clamp(t0, 0., 1.));
        }
        if(t1>=-tol && t1<=1+tol)
        {
          addUnique(times, clamp(t1, 0., 1.));
        }
      }
      
      for (int i=0; i< (int) times.size(); ++i)
      {
        double ti = times[i];
        double val = std::fabs(signed_volume((1-ti)*x0+ti*xnew0,
                                           (1-ti)*x1+ti*xnew1,
                                           (1-ti)*x2+ti*xnew2,
                                           (1-ti)*x3+ti*xnew3));
        errors.push_back(val);
      }
      
      return;
    }
    else
    { // cubic is not monotone: divide up [0,1] accordingly =====================================
      double t0, t1;
      if(B>0)
        t0=(-B-std::sqrt(discriminant))/(3*A);
      else
        t0=(-B+std::sqrt(discriminant))/(3*A);
      t1=C/(3*A*t0);
      if(t1<t0)
        std::swap(t0,t1);
      interval_times.push_back(0);
      if(t0>0 && t0<1)
        interval_times.push_back(t0);
      if(t1>0 && t1<1)
        interval_times.push_back(t1);
      
      interval_times.push_back(1);
    }
  }
  
  // look for roots in indicated intervals ==============================================================
  // evaluate coplanarity more accurately at each endpoint of the intervals
  std::vector<double> interval_values(interval_times.size());
  for(unsigned int i=0; i<interval_times.size(); ++i){
    double t=interval_times[i];
    interval_values[i]=signed_volume((1-t)*x0+t*xnew0, (1-t)*x1+t*xnew1, (1-t)*x2+t*xnew2, (1-t)*x3+t*xnew3);
  }
  // first look for interval endpoints that are close enough to zero, without a sign change
  for(unsigned int i=0; i<interval_times.size(); ++i){
    if(interval_values[i]==0)
    {
      times.push_back(interval_times[i]);
    }
    else if(std::fabs(interval_values[i])<convergence_tol)
    {
      if((i==0 || (interval_values[i-1]>=0 && interval_values[i]>=0) ||
          (interval_values[i-1]<=0 && interval_values[i]<=0)) &&
         (i==interval_times.size()-1 || (interval_values[i+1]>=0 && interval_values[i]>=0) ||
          (interval_values[i+1]<=0 && interval_values[i]<=0)))
      {
        times.push_back(interval_times[i]);
      }
    }
  }
  // and then search in intervals with a sign change
  for(unsigned int i=1; i<interval_times.size(); ++i)
  {
    double tlo=interval_times[i-1], thi=interval_times[i], tmid;
    double vlo=interval_values[i-1], vhi=interval_values[i], vmid;
    if((vlo<0 && vhi>0) || (vlo>0 && vhi<0)){
      // start off with secant approximation (in case the cubic is actually linear)
      double alpha=vhi/(vhi-vlo);
      tmid=alpha*tlo+(1-alpha)*thi;
      for(int iteration=0; iteration<50; ++iteration){
        vmid=signed_volume((1-tmid)*x0+tmid*xnew0, (1-tmid)*x1+tmid*xnew1,
                           (1-tmid)*x2+tmid*xnew2, (1-tmid)*x3+tmid*xnew3);
        if(std::fabs(vmid)<1e-2*convergence_tol) break;
        if((vlo<0 && vmid>0) || (vlo>0 && vmid<0)){ // if sign change between lo and mid
          thi=tmid;
          vhi=vmid;
        }else{ // otherwise sign change between hi and mid
          tlo=tmid;
          vlo=vmid;
        }
        if(iteration%2) alpha=0.5; // sometimes go with bisection to guarantee we make progress
        else alpha=vhi/(vhi-vlo); // other times go with secant to hopefully get there fast
        tmid=alpha*tlo+(1-alpha)*thi;
      }
      times.push_back(tmid);
    }
  }
  std::sort(times.begin(), times.end());
  
  for (int i=0; i< (int) times.size(); ++i)
  {
    double ti = times[i];
    double val = std::fabs(signed_volume((1-ti)*x0+ti*xnew0,
                                       (1-ti)*x1+ti*xnew1,
                                       (1-ti)*x2+ti*xnew2,
                                       (1-ti)*x3+ti*xnew3));
    errors.push_back(val);
  }
}

void BridsonStepper::testCoplanarityTime()
{
  // Inside triangle, time of impact 0.5
  
  // Face positions
  Vec3d xa0(1.0,-1.0,0.0);
  Vec3d xa1(1.0,1.0,0.0);
  Vec3d xa2(-1.0,0.0,0.0);
  // Vert positions
  Vec3d xa3(0.0,0.0,1.0);
  
  // Face velocities
  Vec3d va0(0.0,0.0,0.0);
  Vec3d va1(0.0,0.0,0.0);
  Vec3d va2(0.0,0.0,0.0);
  // Vert velocities
  Vec3d va3(0.0,0.0,-2.0);
  
  // Code takes two sets of vertices, not a set of verts and velocities
  double dta = 50.0;
  
  std::vector<double> timesa; 
  std::vector<double> errorsa;
  getCoplanarityTimes(xa0,xa1,xa2,xa3,xa0+dta*va0,xa1+dta*va1,xa2+dta*va2,xa3+dta*va3,timesa,errorsa);
  
  assert( timesa.size() == 1 );
  assert( fabs(dta*timesa[0]-0.5) < 1.0e-6 );
  
  
  // Outside triangle, time of impact 1.0
  xa3[0] = -20.0; xa3[1] = 100.0; xa3[2] = 2.0;
  getCoplanarityTimes(xa0,xa1,xa2,xa3,xa0+dta*va0,xa1+dta*va1,xa2+dta*va2,xa3+dta*va3,timesa,errorsa);
  
  assert( timesa.size() == 1 );
  assert( fabs(dta*timesa[0]-1.0) < 1.0e-6 );
  
  
  // Other side of triangle. Ensure no roots found
  xa3[0] = 0.0; xa3[1] = 0.0; xa3[2] = -2.0;
  getCoplanarityTimes(xa0,xa1,xa2,xa3,xa0+dta*va0,xa1+dta*va1,xa2+dta*va2,xa3+dta*va3,timesa,errorsa);
  assert( timesa.size() == 0 );
  
  
  // Collision outside of the timestep window.
  double dtbackup = dta;
  dta = 0.1;
  xa3[0] = 10.0; xa3[1] = 1.0; xa3[2] = 50.0;
  va3[0] = -10.0; va3[1] = 15.0; va3[2] = -1.0;
  getCoplanarityTimes(xa0,xa1,xa2,xa3,xa0+dta*va0,xa1+dta*va1,xa2+dta*va2,xa3+dta*va3,timesa,errorsa);
  assert( timesa.size() == 0 );
  dta = dtbackup;
  
  
  // Face moving in orthogonal direction to vert's velocity
  // Vert	
  xa3[0] = 0.0; xa3[1] = 0.0; xa3[2] = 2.0;
  va3[0] = 2.5; va3[1] = 9.0; va3[2] = -0.5;
  // Face
  va0[0] = -1.1; va0[1] = 2.1; va0[2] = 0.0;
  va1[0] = -1.1; va1[1] = 2.1; va1[2] = 0.0;
  va2[0] = -1.1; va2[1] = 2.1; va2[2] = 0.0;
  getCoplanarityTimes(xa0,xa1,xa2,xa3,xa0+dta*va0,xa1+dta*va1,xa2+dta*va2,xa3+dta*va3,timesa,errorsa);
  assert( timesa.size() == 1 );
  assert( fabs(dta*timesa[0]-4.0) < 1.0e-6 );
  
  
  // Vertex not moving
  // Vert	
  xa3[0] = 0.0; xa3[1] = 0.0; xa3[2] = 2.0;
  va3[0] = 0.0; va3[1] = 0.0; va3[2] = 0.0;
  // Face
  va0[0] = -1.1; va0[1] = 2.1; va0[2] = 0.0;
  va1[0] = -1.1; va1[1] = 2.1; va1[2] = 0.0;
  va2[0] = -1.1; va2[1] = 2.1; va2[2] = 0.0;
  getCoplanarityTimes(xa0,xa1,xa2,xa3,xa0+dta*va0,xa1+dta*va1,xa2+dta*va2,xa3+dta*va3,timesa,errorsa);
  assert( timesa.size() == 0 );
  
  
  // Vertex moving away from face
  // Vert	
  xa3[0] = 0.0; xa3[1] = 0.0; xa3[2] = 2.0;
  va3[0] = 0.0; va3[1] = 0.0; va3[2] = 2.0;
  // Face
  va0[0] = -1.1; va0[1] = 2.1; va0[2] = 0.0;
  va1[0] = -1.1; va1[1] = 2.1; va1[2] = 0.0;
  va2[0] = -1.1; va2[1] = 2.1; va2[2] = 0.0;
  getCoplanarityTimes(xa0,xa1,xa2,xa3,xa0+dta*va0,xa1+dta*va1,xa2+dta*va2,xa3+dta*va3,timesa,errorsa);
  assert( timesa.size() == 0 );
  
  
  // Vertex and plane moving towards each other
  // Vert	
  xa3[0] = 0.0; xa3[1] = 0.0; xa3[2] = 2.0;
  va3[0] = 0.0; va3[1] = 0.0; va3[2] = -1.0;
  // Face
  va0[0] = -1.1; va0[1] = 2.1; va0[2] = 1.0;
  va1[0] = -1.1; va1[1] = 2.1; va1[2] = 1.0;
  va2[0] = -1.1; va2[1] = 2.1; va2[2] = 1.0;
  getCoplanarityTimes(xa0,xa1,xa2,xa3,xa0+dta*va0,xa1+dta*va1,xa2+dta*va2,xa3+dta*va3,timesa,errorsa);
  assert( timesa.size() == 1 );
  assert( fabs(dta*timesa[0]-1.0) < 1.0e-6 );
  
  std::cout << "***Vertex Face Coplanarity Tests Successfull" << std::endl;
  
  
  
  
  
  // Edge 0
  Vec3d xb0(-1.0,0.0,0.0);
  Vec3d xb1(1.0,0.0,0.0);
  // Edge 1
  Vec3d xb2(0.0,-1.0,0.0);
  Vec3d xb3(0.0,1.0,0.0);
  
  // Face velocities
  Vec3d vb0(0.0,0.0,0.0);
  Vec3d vb1(0.0,0.0,0.0);
  Vec3d vb2(0.0,0.0,0.0);
  Vec3d vb3(0.0,0.0,0.0);
  
  double dtb = 50.0;
  
  // Never collide - well, kind of. Looks like they are on top of eachother the way these params passed?
  std::vector<double> timesb; 
  std::vector<double> errorsb;
  //getCoplanarityTimes(xb0,xb1,xb2,xb3,xb0+dtb*vb0,xb1+dtb*vb1,xb2+dtb*vb2,xb3+dtb*vb3,timesb,errorsb);
  //for( int i = 0; i < (int) timesb.size(); ++i ) std::cout << timesb[i] << " "; std::cout << std::endl;
  //assert( timesb.size() == 0 );
  
  // Moving apart - never collide
  xb0[0] = -1.0; xb0[1] = 0.0; xb0[2] = 0.0;
  xb1[0] = 1.0; xb1[1] = 0.0; xb1[2] = 0.0;
  xb2[0] = 0.0; xb2[1] = -1.0; xb2[2] = 1.0;
  xb3[0] = 0.0; xb3[1] = 1.0; xb3[2] = 1.0;
  vb0[0] = 0.0; vb0[1] = 0.0; vb0[2] = -1.0;
  vb1[0] = 0.0; vb1[1] = 0.0; vb1[2] = -1.0;
  vb2[0] = 0.0; vb2[1] = 0.0; vb2[2] = 1.0;
  vb3[0] = 0.0; vb3[1] = 0.0; vb3[2] = 1.0;
  getCoplanarityTimes(xb0,xb1,xb2,xb3,xb0+dtb*vb0,xb1+dtb*vb1,xb2+dtb*vb2,xb3+dtb*vb3,timesb,errorsb);
  assert( timesb.size() == 0 );
  
  
  // Edges Moving - collide at time 0.5
  xb0[0] = -1.0; xb0[1] = 0.0; xb0[2] = 0.0;
  xb1[0] = 1.0; xb1[1] = 0.0; xb1[2] = 0.0;
  xb2[0] = 0.0; xb2[1] = -1.0; xb2[2] = 1.0;
  xb3[0] = 0.0; xb3[1] = 1.0; xb3[2] = 1.0;
  vb0[0] = 0.0; vb0[1] = 0.0; vb0[2] = 1.0;
  vb1[0] = 0.0; vb1[1] = 0.0; vb1[2] = 1.0;
  vb2[0] = 0.0; vb2[1] = 0.0; vb2[2] = -1.0;
  vb3[0] = 0.0; vb3[1] = 0.0; vb3[2] = -1.0;
  getCoplanarityTimes(xb0,xb1,xb2,xb3,xb0+dtb*vb0,xb1+dtb*vb1,xb2+dtb*vb2,xb3+dtb*vb3,timesb,errorsb);
  assert( timesb.size() == 1 );
  assert( fabs(dtb*timesb[0]-0.5) < 1.0e-6 );
  
  
  // Edges Moving - collide at time 0.5
  xb0[0] = -1.0; xb0[1] = 0.0; xb0[2] = 0.0;
  xb1[0] = 1.0; xb1[1] = 0.0; xb1[2] = 0.0;
  xb2[0] = 0.0; xb2[1] = -1.0; xb2[2] = 1.0;
  xb3[0] = 0.0; xb3[1] = 1.0; xb3[2] = 1.0;
  vb0[0] = -0.5; vb0[1] = 0.0; vb0[2] = 1.0;
  vb1[0] = 0.5; vb1[1] = 0.0; vb1[2] = 1.0;
  vb2[0] = 0.0; vb2[1] = -0.5; vb2[2] = -1.0;
  vb3[0] = 0.0; vb3[1] = 0.5; vb3[2] = -1.0;
  getCoplanarityTimes(xb0,xb1,xb2,xb3,xb0+dtb*vb0,xb1+dtb*vb1,xb2+dtb*vb2,xb3+dtb*vb3,timesb,errorsb);
  assert( timesb.size() == 1 );
  assert( fabs(dtb*timesb[0]-0.5) < 1.0e-6 );
  
  std::cout << "***Edge Edge Coplanarity Tests Successfull" << std::endl;
}

}





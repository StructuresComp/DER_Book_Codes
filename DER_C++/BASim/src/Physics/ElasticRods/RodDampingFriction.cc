/**
 * \file RodStretchingForce.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/01/2009
 */

#include "BASim/Physics"


using namespace std;

namespace BASim {

RodDampingFriction::RodDampingFriction()
  : m_enbld(false), object_vel(Vec3d(0,0,0)), m_dt(0.01)
  ,m_coulomb_static(0.0), m_coulomb_kinetic(0.0)
  ,m_contact_method(0), m_friction_method(0)
{}

void RodDampingFriction::computeForce(ElasticRod& rod, VecXd& force)
{
	if (m_enbld)
	{
		RodBoundaryCondition* bc = rod.getBoundaryCondition(); 
		int nb = bc->getNumberOfVerticallyScriptedVertices();	// vertices falling under the belt
	  for (int i = 0; i < nb; ++i) {
	  	int vid = bc->getVerticalScriptedVertId(i);	// index in rod class
	  	const Vec3d& vel = rod.getVelocity(vid);
	  	// Normal is always (0, 1, 0)
	  	Vec3d tngt_vel = Vec3d(vel[0], 0.0, vel[2]);
	  	Vec3d rel_vel = tngt_vel - object_vel;
	  	
	  	Scalar normal_force = 0.0;
	  	
			if (m_contact_method == 2)
			{				
				Scalar dv = bc->getDesiredVerticalVelocity(i) - vel[1];
				Scalar df = rod.getVertexMass(vid) * dv / (m_dt);
			
//				Scalar dx = bc->getDesiredVerticalPosition(i) - bc->getUnconstrainedVerticalPosition(i);

//				Scalar x1 = rod.getVertex(vid)[1] + rod.getVelocity(vid)[1] * m_dt;
//				Scalar dx = bc->getDesiredVerticalPosition(i) - x1;
				
//				Scalar dx = bc->getDesiredVerticalPosition(i) - rod.getVertex(vid)[1];
//				Scalar df = rod.getVertexMass(vid) * (dx / (m_dt * m_dt));

				force(rod.vertIdx(vid, 1)) += df;			
				normal_force = df;
			} else {
				Scalar unconsPos = bc->getUnconstrainedVerticalPosition(i);
				Scalar dy = bc->getDesiredVerticalPosition(i) - unconsPos;
				normal_force = rod.getVertexMass(vid) * dy / (m_dt * m_dt);
			}
				  	
	  	Vec3d slide_friction(0,0,0); // = - m_friction * rel_vel;	// damping friction, -cv
	  	
	  	if (m_friction_method == 5) {
	  		Vec3d friction_force = Vec3d(0,0,0);
	  		bool is_static = bc->isStatic(i);
	  		Scalar coulomb_limit = 0.0;

				Scalar unconsPos = bc->getUnconstrainedVerticalPosition(i);
				Scalar dy = bc->getDesiredVerticalPosition(i) - unconsPos;
				Scalar normalF = rod.getVertexMass(vid) * dy / (m_dt * m_dt);	 
				
				// FOR TEST
				// This is a velocity, at the end of current time step, which includes the other forces.
		  	Vec3d rel_vel1 = rel_vel;
		  	Vec3d other_force = force.segment(rod.vertIdx(vid, 0), 3);
		  	
		  	other_force[1] = 0;
		  	rel_vel1 += other_force / rod.getVertexMass(vid) * m_dt;
				
//	  		std::cout << normalF << " " << is_static << " -> ";
	  		
	  		if (is_static || rel_vel.norm() < 1e-6) {
					coulomb_limit = m_coulomb_static * fabs(normalF);
	  			is_static = true;

	  			// Kill all relative tangential velocity
//	  			friction_force = - rod.getVertexMass(vid) * rel_vel1 / m_dt;
	  			friction_force = - rod.getVertexMass(vid) * rel_vel / m_dt;
	  			
	  			if (friction_force.norm() > coulomb_limit) {
	  				is_static = false;
	  			}
	  		} 
	  		
	  		if (!is_static) {
					coulomb_limit = m_coulomb_kinetic * fabs(normalF);
	  			if (coulomb_limit >= rod.getVertexMass(vid) * rel_vel.norm() / m_dt) {
	  				is_static = true;
	  				friction_force = - rod.getVertexMass(vid) * rel_vel / m_dt;
	  			} else {
		  			friction_force = - coulomb_limit * rel_vel / rel_vel.norm();
		  		}
	  		}
	  		
//	  		std::cout << is_static << " " << friction_force << " " << force.segment(rod.vertIdx(vid, 0), 3) << rel_vel <<  " \n";
	  		
	  		bc->setStatic(i, is_static);
				slide_friction = friction_force;
	  	}
	  	
	  	// CASE 6 :: Always static friction + Coulomb's limit (with static coefficient)
	  	if (m_friction_method == 6) {
	  		Vec3d friction_force = Vec3d(0,0,0);
	  		Scalar coulomb_limit = 0.0;

				//normal_force				
//	  		std::cout << normalF << " " << is_static << " -> ";
	  		
				coulomb_limit = m_coulomb_static * fabs(normal_force);

  			// Kill all relative tangential velocity
  			friction_force = - rod.getVertexMass(vid) * rel_vel / m_dt;
  			
  			if (friction_force.norm() > coulomb_limit) {
  				if (rel_vel.norm() > 1e-6) {
	  				friction_force = - coulomb_limit * rel_vel / rel_vel.norm();
	  			} else {
	  				friction_force = Vec3d(0,0,0);
	  			}
  			}

				slide_friction = friction_force;
	  	}

	  	if (m_friction_method == 7) {
	  		Vec3d friction_force = Vec3d(0,0,0);
	  		Scalar coulomb_limit = m_coulomb_static * fabs(normal_force);

  			// Kill all relative tangential velocity
  			friction_force = - rod.getVertexMass(vid) * rel_vel / m_dt;
  			
	  		bc->setStatic(i, false);
	  		  			
  			if (friction_force.norm() > coulomb_limit) {
  				if (rel_vel.norm() > 1e-6) {
	  				friction_force = - coulomb_limit * rel_vel / rel_vel.norm();
			  		bc->setStatic(i, true);
	  			} else {
	  				friction_force = Vec3d(0,0,0);
	  			}
  			}

				slide_friction = friction_force;
	  	}
	  		  	
	  	force.segment(rod.vertIdx(vid, 0), 3) += slide_friction;
	  }
	}
}

void RodDampingFriction::computeForceDX(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J)
{
  return;
  
	if (m_enbld)
	{
		if (m_contact_method == 2)
		{	
			const RodBoundaryCondition* bc = rod.getBoundaryCondition(); 
			int nb = bc->getNumberOfVerticallyScriptedVertices();
			for (int i = 0; i < nb; ++i) {
				int vid = bc->getVerticalScriptedVertId(i);

				{
					Scalar df = - rod.getVertexMass(vid) / (m_dt * m_dt);
			
//					J.add(baseidx + rod.vertIdx(vid, 1), baseidx + rod.vertIdx(vid, 1), df * scale);
				}
			}
		}
	}
}

void RodDampingFriction::computeForceDV(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J)
{
	if (m_enbld)
	{
		const RodBoundaryCondition* bc = rod.getBoundaryCondition(); 
		int nb = bc->getNumberOfVerticallyScriptedVertices();
	  for (int i = 0; i < nb; ++i) {
	  	int vid = bc->getVerticalScriptedVertId(i);

			if (m_contact_method == 2) {
				Scalar df = - rod.getVertexMass(vid) / (m_dt);
				J.add(baseidx + rod.vertIdx(vid, 1), baseidx + rod.vertIdx(vid, 1), df * scale);
			}

      // Friction			      
      if (m_friction_method != 6 && m_friction_method != 7) continue;
      
//	  	Scalar fdv = - m_friction;

//			continue;
			
	  	// TODO: Set upper-bound, 
	 	
	  	if (m_friction_method == 5) {
	  		bool is_static = bc->isStatic(i);

	  		if (is_static) {
	  			Scalar df = - rod.getVertexMass(vid) / m_dt;
	  			
				  J.add(rod.vertIdx(vid, 0), rod.vertIdx(vid, 0), df * scale);
				  J.add(rod.vertIdx(vid, 2), rod.vertIdx(vid, 2), df * scale);
	  		} else {
					const Vec3d& vel = rod.getVelocity(vid);
					Vec3d tngt_vel = Vec3d(vel[0], 0.0, vel[2]);
					Vec3d rel_vel = tngt_vel - object_vel;
	  		  		
					Scalar unconsPos = bc->getUnconstrainedVerticalPosition(i);
					Scalar dy = bc->getDesiredVerticalPosition(i) - unconsPos;
					Scalar normalF = fabs (rod.getVertexMass(vid) * dy / (m_dt * m_dt));	 
					
					Scalar vnorm = rel_vel.norm();
					Scalar vnorm3 = vnorm * vnorm * vnorm;
					
				  J.add(rod.vertIdx(vid, 0), rod.vertIdx(vid, 0), - m_coulomb_kinetic * normalF * scale * ( 1.0 / vnorm - rel_vel[0] * rel_vel[0] / vnorm3) );
				  J.add(rod.vertIdx(vid, 2), rod.vertIdx(vid, 2), - m_coulomb_kinetic * normalF * scale * ( 1.0 / vnorm - rel_vel[2] * rel_vel[2] / vnorm3) );
				  J.add(rod.vertIdx(vid, 0), rod.vertIdx(vid, 2), - m_coulomb_kinetic * normalF * scale * (             - rel_vel[0] * rel_vel[2] / vnorm3) );
				  J.add(rod.vertIdx(vid, 2), rod.vertIdx(vid, 0), - m_coulomb_kinetic * normalF * scale * (             - rel_vel[0] * rel_vel[2] / vnorm3) );
	
	  		}
				
	  	}
	  	else if (m_friction_method == 6) {
  			Scalar df = - rod.getVertexMass(vid) / m_dt;
  			
			  J.add(rod.vertIdx(vid, 0), rod.vertIdx(vid, 0), df * scale);
			  J.add(rod.vertIdx(vid, 2), rod.vertIdx(vid, 2), df * scale);
	  	}
	  	else if (m_friction_method == 7) {
	  	
//		    J.add(rod.vertIdx(vid, 0), rod.vertIdx(vid, 0), -0.3 * scale);
//		    J.add(rod.vertIdx(vid, 2), rod.vertIdx(vid, 2), -0.3 * scale);

		    J.add(rod.vertIdx(vid, 0), rod.vertIdx(vid, 0), -m_coulomb_static * scale);
		    J.add(rod.vertIdx(vid, 2), rod.vertIdx(vid, 2), -m_coulomb_static * scale);
	  	}



	  }
	}
}


} // namespace BASim

/**
 * \file RodDampingFriction.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/08/2009
 */

#ifndef RODFRICTIONDAMPING_HH
#define RODFRICTIONDAMPING_HH

namespace BASim {

/** This class implements a mass damping force for a rod. */
class RodDampingFriction : public RodExternalForce
{
public:

  /** Constructor for creating the mass damping force

      \param[in] damping The damping coefficient.
  */
  RodDampingFriction();
  
  void setObjectVel(const Vec3d& v) {
  	object_vel = v;
  	object_vel[1] = 0.0;
  }
  
  void setTimeStep( Scalar dt ) {
  	m_dt = dt;
  }
  
  void enable()
  {
  	m_enbld = true;
  }
  
  void disable()
  {
  	m_enbld = false;
  }

  void setContactAndFrictionMethod( int contact, int friction ) 
  {
		m_contact_method = contact;
		m_friction_method = friction;    	
  }
    
  /** Computes the mass damping force for the given rod and adds it to
      the given vector.

      \param[in] rod The rod to compute the force for.
      \param[out] force Vector storing the forces on the rod.
  */
  virtual void computeForce(ElasticRod& rod, VecXd& force);
//  virtual void computeForce(const ElasticRod& rod, VecXd& force);

  /** Computes the derivative of the mass damping force with respect
      to the velocities of the rod. */
  virtual void computeForceDV(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J);

  virtual void computeForceDX(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J);

  void setCoulombCoefficient( double static_cof, double kinetic_cof ) {
  	m_coulomb_static = static_cof;
  	m_coulomb_kinetic = kinetic_cof;
  }

protected:

  bool m_enbld;
  Vec3d object_vel;	// only using x and z components
  Scalar m_dt;

  Scalar m_coulomb_static;
  Scalar m_coulomb_kinetic;
  
  
	int m_contact_method;
	int m_friction_method;    
  
};

} // namespace BASim

#endif // RodDampingFriction_HH

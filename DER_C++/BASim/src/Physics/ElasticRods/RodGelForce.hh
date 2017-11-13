/**
 * \file RodGelForce.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 05/13/2010
 */

#ifndef RODGELFORCE_HH
#define RODGELFORCE_HH

#include "BASim/src/Physics/ElasticRods/RodExternalForce.hh"

namespace BASim {

/** External force that simulates the resistance of an elastic medium
    surrounding the rod. */
class RodGelForce : public RodExternalForce
{
public:

  /**
   * Constructor for creating the gel force.
   *
   * \param[in] stiffness The stiffness of the gel force
   */

  RodGelForce(Scalar stiffness) : m_stiffness(stiffness) {}

  /**
   * Computes the force for the given rod and adds it to the given
   * vector.
   *
   * \param[in] rod The rod to compute the gravity force for.
   * \param[out] force Vector storing the forces on the rod.
   */

//  void computeForce(const ElasticRod& rod, VecXd& force)
  void computeForce(ElasticRod& rod, VecXd& force)
  {
    // need implementation of the force
  }

  void computeForceDX(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J)
  {
    // need implementation of the Jacobian of the force with respect
    // to positions
  }

  void computeForceDV(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J)
  {
    // need implementation of the Jacobian of the force with respect
    // to velocities
  }

protected:

  Scalar m_stiffness;
};

} // namespace BASim

#endif // RODGELFORCE_HH

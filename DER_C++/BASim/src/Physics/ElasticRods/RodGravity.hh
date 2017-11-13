/**
 * \file RodGravity.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/04/2009
 */

#ifndef RODGRAVITY_HH
#define RODGRAVITY_HH

#include "BASim/src/Physics/ElasticRods/RodExternalForce.hh"

namespace BASim {

/** External force that applies gravity to a rod. */
class RodGravity : public RodExternalForce
{
public:

  /**
   * Constructor for creating the gravity force.
   *
   * \param[in] gravity The acceleration due to gravity to be applied to rods.
   */

  RodGravity(const Vec3d& gravity)
    : m_gravity(gravity)
  {}

  /**
   * Computes the force due to gravity for the given rod and adds it
   * to the given vector.
   *
   * \param[in] rod The rod to compute the gravity force for.
   * \param[out] force Vector storing the forces on the rod.
   */

//  void computeForce(const ElasticRod& rod, VecXd& force)
  void computeForce(ElasticRod& rod, VecXd& force)
  {
    for (int i = 0; i < rod.nv(); ++i) {
      Vec3d f = rod.getVertexMass(i) * m_gravity;
      for (int coord = 0; coord < 3; ++coord) {
        force(rod.vertIdx(i, coord)) += f(coord);
      }
    }
  }

  void computeForceDX(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J) {}

  void computeForceDV(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J) {}


protected:

  Vec3d m_gravity;
};

} // namespace BASim

#endif // RODGRAVITY_HH

/**
 * \file RodMassDamping.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/08/2009
 */

#ifndef RODMASSDAMPING_HH
#define RODMASSDAMPING_HH

#include "BASim/src/Physics/ElasticRods/RodExternalForce.hh"

namespace BASim {

/** This class implements a mass damping force for a rod. */
class RodMassDamping : public RodExternalForce
{
public:

  /** Constructor for creating the mass damping force

      \param[in] damping The damping coefficient.
  */
  RodMassDamping(const Scalar& damping)
    : m_damping(damping)
  {}

  /** Computes the mass damping force for the given rod and adds it to
      the given vector.

      \param[in] rod The rod to compute the force for.
      \param[out] force Vector storing the forces on the rod.
  */
//  void computeForce(const ElasticRod& rod, VecXd& force)
  void computeForce( ElasticRod& rod, VecXd& force)
  {
    for (int i = 0; i < rod.ndof(); ++i) {
      force(i) -= m_damping * rod.getMass(i) * rod.getVel(i);
    }
  }

  void computeForceDX(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J) {}

  /** Computes the derivative of the mass damping force with respect
      to the velocities of the rod. */
  void computeForceDV(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J)
  {
    for (int i = 0; i < rod.ndof(); ++i) {
      J.add(baseidx+i, baseidx+i, -m_damping * rod.getMass(i) * scale);
    }
  }

protected:

  Scalar m_damping;
};

} // namespace BASim

#endif // RODMASSDAMPING_HH

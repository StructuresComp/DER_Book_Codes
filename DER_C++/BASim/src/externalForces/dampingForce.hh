/**
 * \file dampingForce.hh
 *
 * \author Khalid Jawed
 *
 * Built upon RodMassDamping.hh
 */

#ifndef DAMPINGFORCE_HH
#define DAMPINGFORCE_HH

#include "BASim/src/Physics/ElasticRods/RodExternalForce.hh"

namespace BASim
{

/** This class implements a damping force for flagella. */
class dampingForce: public RodExternalForce
{
public:

    /** Constructor for creating the damping force

     \param[in] damping The damping coefficient.
     */
    dampingForce(const Scalar& eta, Scalar & m_time) :
        m_eta(eta)
    {
        ttime = &m_time;
    }


    void computeForce(ElasticRod& rod, VecXd& force)
    {
        for (int i = 0; i < rod.nv(); ++i)
        {
            Vec3d t = rod.getTangent(i);
            Vec3d v = rod.getVelocity(i);
            Vec3d f = - m_eta * v;
            force.segment<3>(i * 4) += f * rod.getVoronoiLength(i);
        }
        
    }

    void computeForceDX(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J)
    {
		;
    }

    /** Computes the derivative of the mass damping force with respect
     to the velocities of the rod. */
    void computeForceDV(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J)
    {        
        for (int i = 0; i < rod.nv(); ++i)
        {
            Vec3d t = rod.getTangent(i);
            Vec3d v = rod.getVelocity(i);
            Mat3d j = - m_eta * Mat3d::Identity();
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    J.add(baseidx + i * 4 + k, baseidx + i * 4 + l, scale * j(k, l) * rod.getVoronoiLength(i));
        }
    }

protected:

    Scalar m_eta;
    Scalar *ttime;

};

} // namespace BASim

#endif // FLAGELLADAMPING_HH

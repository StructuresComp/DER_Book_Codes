/**
 * \file terminalLoads.hh
 *
 * \author Khalid Jawed
 * \date 06/25/2017
 *
 */

#ifndef TERMINALLOADS_HH
#define TERMINALLOADS_HH
#include "BASim/src/Physics/ElasticRods/RodExternalForce.hh"

namespace BASim
{
    class terminalLoads : public RodExternalForce
    {
    public:
        
        /** Constructor
         */
        terminalLoads(Scalar k1i, Scalar k2i, Scalar t1i, 
          Scalar k1f, Scalar k2f, Scalar t1f,  
          Scalar EI, Scalar GJ, Scalar & m_time)
        {
			v1_i = k1i;
			v2_i = k2i;
			v3_i = t1i;
			v1_f = k1f;
			v2_f = k2f;
			v3_f = t1f;
			m_EI = EI;
			m_GJ = GJ;
            ttime = &m_time;
        }        
        
        void computeForce(ElasticRod& rod, VecXd& force)
        {
            int N = rod.nv();
            int nDof = 4*N - 1;
            
            Scalar multiplier1 = *ttime / 10.0;
            if (multiplier1 > 1.0) multiplier1 = 1.0;
            
            Scalar multiplier2 = *ttime / 10.0;
            if (multiplier2 > 1.0) multiplier2 = 1.0;

            Vec3d m1i = rod.getMaterial1(0);
            Vec3d m2i = rod.getMaterial2(0);
            Vec3d m1f = rod.getMaterial1(N-2);
            Vec3d m2f = rod.getMaterial2(N-2);
            Scalar dli = rod.getEdgeLength(0);
            Scalar dlf = rod.getEdgeLength(N-2);
            
            // Torsion
            Scalar mT_i = - m_GJ * v3_i; // i=0
            Scalar mT_f =   m_GJ * v3_f; // i=N-2
            force( 3 ) += mT_i * multiplier1;
            force( nDof - 4 ) += mT_f * multiplier1;
            
            // Forces
            // Node 0
            Vec3d f0 = m_EI / dli * ( - v1_i * m2i + v2_i * m1i );
            // Node 1
            Vec3d f1 = - f0;
            // Node N-2
            Vec3d fN1 = m_EI / dlf * ( v1_f * m2f - v2_f * m1f );
            // Node N-1
            Vec3d fN = - fN1;
            
            force.segment<3>(0 * 4) += f0 * multiplier2;
            force.segment<3>(1 * 4) += f1 * multiplier2;
            force.segment<3>((N-2) * 4) += fN1 * multiplier2;
            force.segment<3>((N-1) * 4) += fN * multiplier2;
        }
        
        void computeForceDX(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J) {}
        
        void computeForceDV(int baseidx, const ElasticRod& rod, Scalar scale, MatrixBase& J) {}

    protected:
        Scalar v1_i, v2_i, v3_i, v1_f, v2_f, v3_f;
        Scalar m_EI, m_GJ;
        Scalar *ttime;
    };
    
} // namespace BASim

#endif // STOKESLET_HH


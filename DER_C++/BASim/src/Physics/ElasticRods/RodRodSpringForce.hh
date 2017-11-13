/**
 * \file RodRodSpringForce.hh
 *
 * \author smith@cs.columbia.edu
 * \date 06/02/2010
 */

#ifndef RODRODSPRINGFORCE_HH
#define RODRODSPRINGFORCE_HH

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Physics/ElasticRods/RodRodExternalForce.hh"

namespace BASim
{

class RodRodSpringForce : public RodRodExternalForce
{
public:

  RodRodSpringForce( const std::vector<ElasticRod*>& rods, int rodA, int rodB, int vertA, int vertB, double k, double l0 )
  : RodRodExternalForce(true)
  , m_rods(rods)
  , m_rod_A(rodA)
  , m_rod_B(rodB)
  , m_vrt_A(vertA)
  , m_vrt_B(vertB)
  , m_k(k)
  , m_l0(l0)
  , m_base_idx_A(0)
  , m_base_idx_B(0)
  {
    assert( m_rod_A >= 0 );
    assert( m_rod_A < (int) m_rods.size() );
    assert( m_rod_B >= 0 );
    assert( m_rod_B < (int) m_rods.size() );
    assert( k >= 0.0 );
    assert( l0 >= 0.0 );
    
    for( int i = 0; i < m_rod_A; ++i ) m_base_idx_A += m_rods[i]->ndof();
    for( int i = 0; i < m_rod_B; ++i ) m_base_idx_B += m_rods[i]->ndof();
  }

  virtual ~RodRodSpringForce() {}

  virtual void computeForce( VecXd& force )
  {
    Vec3d x1 = m_rods[m_rod_A]->getVertex(m_vrt_A);
    Vec3d x0 = m_rods[m_rod_B]->getVertex(m_vrt_B);
    
    Vec3d nhat = x1 - x0;
    double l = nhat.norm();
    assert( l != 0.0 );
    nhat /= l;
    
    Vec3d f = -m_k*(l-m_l0)*nhat;
    
    for( int coord = 0; coord < 3; ++coord ) 
    {
      force(m_base_idx_A + m_rods[m_rod_A]->vertIdx(m_vrt_A,coord)) += f(coord);
    }

    for( int coord = 0; coord < 3; ++coord ) 
    {
      force(m_base_idx_B + m_rods[m_rod_B]->vertIdx(m_vrt_B,coord)) += -f(coord);
    }
    
    //std::cout << f.transpose() << std::endl;
  }

  virtual void computeForceDX( Scalar scale, MatrixBase& J )
  {
    Vec3d x1 = m_rods[m_rod_A]->getVertex(m_vrt_A);
    Vec3d x0 = m_rods[m_rod_B]->getVertex(m_vrt_B);

    Vec3d nhat = x1 - x0;
    double l = nhat.norm();
    assert( l != 0.0 );
    nhat /= l;

    Mat3d jac = m_k*nhat*nhat.transpose() + m_k*(l-m_l0)*(Mat3d::Identity()-nhat*nhat.transpose())/l;
    jac *= -scale;

    // x1,x1
    std::vector<int> indicesA;
    for( int coord = 0; coord < 3; ++coord ) 
    {
      indicesA.push_back(m_base_idx_A+m_rods[m_rod_A]->vertIdx(m_vrt_A,coord));
    }
    J.add(indicesA,indicesA,jac);

    // x0,x0
    std::vector<int> indicesB;
    for( int coord = 0; coord < 3; ++coord ) 
    {
      indicesB.push_back(m_base_idx_B+m_rods[m_rod_B]->vertIdx(m_vrt_B,coord));
    }
    J.add(indicesB,indicesB,jac);
    
    // x1,x0
    J.add(indicesA,indicesB,jac);
    // x0,x1
    J.add(indicesB,indicesA,jac);
  }

  virtual void computeForceDV( Scalar scale, MatrixBase& J )
  {
    // No v dependence
  }

private:

  const std::vector<ElasticRod*>& m_rods;
  int m_rod_A;
  int m_rod_B;
  int m_vrt_A;
  int m_vrt_B;
  int m_k;
  int m_l0;
  
  int m_base_idx_A;
  int m_base_idx_B;
};

} // namespace BASim

#endif // RODRODSPRINGFORCE_HH

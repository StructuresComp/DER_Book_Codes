/**
 * \file MinimalRodStateBackup.cc
 *
 * \author smith@cs.columbia.edu
 * \date 07/12/2010
 */

#include "MinimalRodStateBackup.hh"

namespace BASim 
{

void MinimalRodStateBackup::resize( const ElasticRod& rod )
{
  int nv = rod.nv();
  int ne = rod.ne();

  m_vertexPositions.resize(nv);
  m_vertexVelocities.resize(nv);
  m_referenceTwist.resize(nv);

  m_theta.resize(ne);
  m_thetaDot.resize(ne);
  m_referenceDirectors1.resize(ne);
}

void MinimalRodStateBackup::backupRod( ElasticRod& rod )
{
  assert( (int) m_vertexPositions.size() == rod.nv() );
  assert( (int) m_vertexVelocities.size() == rod.nv() );
  assert( (int) m_referenceTwist.size() == rod.nv() );

  assert( (int) m_theta.size() == rod.ne() );
  assert( (int) m_thetaDot.size() == rod.ne() );
  assert( (int) m_referenceDirectors1.size() == rod.ne() );

  // Backup quantities associated with vertices
  for( ElasticRod::vertex_iter itr = rod.vertices_begin(); itr != rod.vertices_end(); ++itr )
    m_vertexPositions[(*itr).idx()] = rod.getVertex(*itr);

  for( ElasticRod::vertex_iter itr = rod.vertices_begin(); itr != rod.vertices_end(); ++itr )
    m_vertexVelocities[(*itr).idx()] = rod.getVelocity((*itr).idx());

  for( ElasticRod::vertex_iter itr = rod.vertices_begin(); itr != rod.vertices_end(); ++itr )
    m_referenceTwist[(*itr).idx()] = rod.getReferenceTwist(*itr);

  // Backup quantities associated with edges
  for( ElasticRod::edge_iter itr = rod.edges_begin(); itr != rod.edges_end(); ++itr )
    m_theta[(*itr).idx()] = rod.getTheta(*itr);

  for( ElasticRod::edge_iter itr = rod.edges_begin(); itr != rod.edges_end(); ++itr )
    m_thetaDot[(*itr).idx()] = rod.getThetaDot((*itr).idx());
  
  for( ElasticRod::edge_iter itr = rod.edges_begin(); itr != rod.edges_end(); ++itr )
    m_referenceDirectors1[(*itr).idx()] = rod.getReferenceDirector1(*itr);
}

void MinimalRodStateBackup::restoreRod( ElasticRod& rod )
{
  assert( (int) m_vertexPositions.size() == rod.nv() );
  assert( (int) m_vertexVelocities.size() == rod.nv() );
  assert( (int) m_referenceTwist.size() == rod.nv() );

  assert( (int) m_theta.size() == rod.ne() );
  assert( (int) m_thetaDot.size() == rod.ne() );
  assert( (int) m_referenceDirectors1.size() == rod.ne() );

  // Restore quantities associated with vertices
  for( ElasticRod::vertex_iter itr = rod.vertices_begin(); itr != rod.vertices_end(); ++itr )
    rod.setVertex(*itr, m_vertexPositions[(*itr).idx()]);

  for( ElasticRod::vertex_iter itr = rod.vertices_begin(); itr != rod.vertices_end(); ++itr )
    rod.setVelocity((*itr).idx(), m_vertexVelocities[(*itr).idx()]);
  
  for( ElasticRod::vertex_iter itr = rod.vertices_begin(); itr != rod.vertices_end(); ++itr )
    rod.setReferenceTwist(*itr, m_referenceTwist[(*itr).idx()]);
  
  // Restore quantities associated with edges
  for( ElasticRod::edge_iter itr = rod.edges_begin(); itr != rod.edges_end(); ++itr )
    rod.setTheta(*itr, m_theta[(*itr).idx()]);

  for( ElasticRod::edge_iter itr = rod.edges_begin(); itr != rod.edges_end(); ++itr )
    rod.setThetaDot((*itr).idx(), m_thetaDot[(*itr).idx()]);
  
  for( ElasticRod::edge_iter itr = rod.edges_begin(); itr != rod.edges_end(); ++itr )
    rod.setReferenceDirector1(*itr, m_referenceDirectors1[(*itr).idx()]);

  rod.computeEdges();
  rod.computeTangents();
  rod.computeCurvatureBinormals();
  rod.computeEdgeLengths();
  rod.computeVoronoiLengths();

  for( ElasticRod::edge_iter itr = rod.edges_begin(); itr != rod.edges_end(); ++itr )
    rod.setReferenceDirector2(*itr, rod.getTangent(*itr).cross(rod.getReferenceDirector1(*itr)));
  
  rod.computeMaterialDirectors();

  rod.updateForceProperties();

  // NOTE: Gradient and hessian should be invalid, but thats ok because the valid flags should be false
  
  // NOTE: Rod undeformed strain should be inavlid, but the first call to viscousUpdate() in the implicit
  //       solver will fix this, as long as the corresponding deformed properties are fine.
  //    viscous stretching ref length
  //    viscous undeformed material curvature vector
  //    viscous undeformed twist
  // BIGBIGBIG: If we add in a line search or something, we might have to back these guys up. Right now
  //            the backup/restore only happens outside of a solve.
  
  #ifdef DEBUG
    ObjPropHandle<bool> vh;
    rod.property_handle(vh,"grad twist valid");
    assert( !rod.property(vh) );
    rod.property_handle(vh,"hess twist valid");
    assert( !rod.property(vh) );
    rod.property_handle(vh,"grad kappa valid");
    assert( !rod.property(vh) );
    rod.property_handle(vh,"hess kappa valid");
    assert( !rod.property(vh) );
    rod.verifyProperties();
  #endif
}

void MinimalRodStateBackup::clear()
{
  m_vertexPositions.clear();
  m_vertexVelocities.clear();
  m_referenceTwist.clear();

  m_theta.clear();
  m_thetaDot.clear();
  m_referenceDirectors1.clear();
}

} // namespace BASim



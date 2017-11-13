/**
 * \file MinimalRodStateBackup.hh
 *
 * \author smith@cs.columbia.edu
 * \date 07/12/2010
 */


#ifndef MINIMALRODSTATEBACKUP_HH
#define MINIMALRODSTATEBACKUP_HH

#include <fstream>

#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Core/Property.hh"

namespace BASim 
{

class MinimalRodStateBackup
{
public:

  void resize( const ElasticRod& rod );

  void backupRod( ElasticRod& rod );

  void restoreRod( ElasticRod& rod );

  void clear();

  // Quantities associated with vertices
  std::vector<Vec3d> m_vertexPositions;
  std::vector<Vec3d> m_vertexVelocities;
  std::vector<Scalar> m_referenceTwist;

  // Quantities associated with edges
  std::vector<Scalar> m_theta;
  std::vector<Scalar> m_thetaDot;
  std::vector<Vec3d> m_referenceDirectors1;
};

} // namespace BASim

#endif // MINIMALRODSTATEBACKUP_HH


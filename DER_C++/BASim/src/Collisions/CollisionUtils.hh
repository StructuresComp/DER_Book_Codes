/**
 * \file CollisionUtils.hh
 *
 * \author (The Internets + Books)
 * \date 04/17/2010
 */

#ifndef COLLISIONUTILS_HH
#define COLLISIONUTILS_HH

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Math/Math.hh"

namespace BASim {

// Closest point on a triangle to a vertex
Vec3d ClosestPtPointTriangle( const Vec3d& p, const Vec3d& a, const Vec3d& b, const Vec3d& c );

// Computes the distance between and closest points of two edges. 
double ClosestPtSegmentSegment( const Vec3d& p1, const Vec3d& q1, const Vec3d& p2, const Vec3d& q2, double& s, double& t, Vec3d& c1, Vec3d& c2 );
  
// Computes the barycentric coordiantes of a point wrt a triangle
void Barycentric( const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& p, double& u, double& v, double& w );

}
  
#endif

// TODO:
//   o It would be nice to handle degenerate cases better in these methods.  
//     They all handle degenerate cases, but getting PREDICTABLE behavior out would rock!!!

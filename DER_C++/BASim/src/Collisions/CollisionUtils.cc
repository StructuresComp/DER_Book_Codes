#include "CollisionUtils.hh"

namespace BASim 
{

// Adapted from Christer Ericson, "Real Time Collision Detection"
Vec3d ClosestPtPointTriangle( const Vec3d& p, const Vec3d& a, const Vec3d& b, const Vec3d& c )
{
  // Check if P in vertex region outside A
  const Vec3d ab = b - a;
  const Vec3d ac = c - a;
  const Vec3d ap = p - a;
  double d1 = ab.dot(ap);
  double d2 = ac.dot(ap);
  if (d1 <= 0.0 && d2 <= 0.0) return a; // barycentric coordinates (1,0,0)
  
  // Check if P in vertex region outside B
  const Vec3d bp = p - b;
  double d3 = ab.dot(bp);
  double d4 = ac.dot(bp);
  if (d3 >= 0.0 && d4 <= d3) return b; // barycentric coordinates (0,1,0)
  
  // Check if P in edge region of AB, if so return projection of P onto AB
  double vc = d1*d4 - d3*d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    double v = d1 / (d1 - d3);
    return a + v * ab; // barycentric coordinates (1-v,v,0)
  }
  
  // Check if P in vertex region outside C
  const Vec3d cp = p - c;
  double d5 = ab.dot(cp);
  double d6 = ac.dot(cp);
  if (d6 >= 0.0 && d5 <= d6) return c; // barycentric coordinates (0,0,1)
  
  // Check if P in edge region of AC, if so return projection of P onto AC
  double vb = d5*d2 - d1*d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    double w = d2 / (d2 - d6);
    return a + w * ac; // barycentric coordinates (1-w,0,w)
  }
  
  // Check if P in edge region of BC, if so return projection of P onto BC
  double va = d3*d6 - d5*d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return b + w * (c - b); // barycentric coordinates (0,1-w,w)
  }
  
  // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
  double denom = 1.0 / (va + vb + vc);
  double v = vb * denom;
  double w = vc * denom;
  return a + ab * v + ac * w; // = u*a + v*b + w*c, u = va * denom = 1.0f - v - w
}

  
// Adapted from Christer Ericson, "Real Time Collision Detection"
// Computes closest points C1 and C2 of S1(s)=P1+s*(Q1-P1) and
// S2(t)=P2+t*(Q2-P2), returning s and t. Function result is squared
// distance between between S1(s) and S2(t).
// TODO: Explore behavior in degenerate case more closely.
double ClosestPtSegmentSegment( const Vec3d& p1, const Vec3d& q1, const Vec3d& p2, const Vec3d& q2, double& s, double& t, Vec3d& c1, Vec3d& c2 )
{
  double EPSILON = 1.0e-9;
  
  Vec3d d1 = q1 - p1; // Direction vector of segment S1
  Vec3d d2 = q2 - p2; // Direction vector of segment S2
  Vec3d r = p1 - p2;
  double a = d1.dot(d1); // Squared length of segment S1, always nonnegative
  double e = d2.dot(d2); // Squared length of segment S2, always nonnegative
  double f = d2.dot(r);
  
  // Check if either or both segments degenerate into points
  if (a <= EPSILON && e <= EPSILON) 
  {
    // Both segments degenerate into points
    s = t = 0.0;
    c1 = p1;
    c2 = p2;
    return (c1-c2).dot(c1-c2);
  }
  if (a <= EPSILON) 
  {
    // First segment degenerates into a point
    s = 0.0;
    t = f / e; // s = 0 => t = (b*s + f) / e = f / e
    t = clamp(t, 0.0, 1.0);
  }
  else 
  {
    double c = d1.dot(r);
    if (e <= EPSILON) 
    {
      // Second segment degenerates into a point
      t = 0.0;
      s = clamp(-c/a, 0.0, 1.0); // t = 0 => s = (b*t - c) / a = -c / a
    } 
    else 
    {
      // The general nondegenerate case starts here
      double b = d1.dot(d2);
      double denom = a*e-b*b; // Always nonnegative
      
      // If segments not parallel, compute closest point on L1 to L2, and
      // clamp to segment S1. Else pick arbitrary s (here 0)
      if (denom != 0.0) 
      {
        s = clamp((b*f - c*e) / denom, 0.0, 1.0);
      } 
      else s = 0.0;
      
      // Compute point on L2 closest to S1(s) using
      // t = Dot((P1+D1*s)-P2,D2) / Dot(D2,D2) = (b*s + f) / e
      t = (b*s + f) / e;
      
      // If t in [0,1] done. Else clamp t, recompute s for the new value
      // of t using s = Dot((P2+D2*t)-P1,D1) / Dot(D1,D1)= (t*b - c) / a
      // and clamp s to [0, 1]
      if (t < 0.0) 
      {
        t = 0.0;
        s = clamp(-c / a, 0.0, 1.0);
      } 
      else if (t > 1.0) 
      {
        t = 1.0;
        s = clamp((b - c) / a, 0.0, 1.0);
      }
    }
  }
  
  c1 = p1 + d1 * s;
  c2 = p2 + d2 * t;
  return (c1 - c2).dot(c1 - c2);
}

// Adapted from Christer Ericson, "Real Time Collision Detection"
// Compute barycentric coordinates (u, v, w) for 
// point p with respect to triangle (a, b, c)
void Barycentric( const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& p, double& u, double& v, double& w )
{
  const Vec3d v0 = b - a;
  const Vec3d v1 = c - a;
  const Vec3d v2 = p - a;
  double d00 = v0.dot(v0);
  double d01 = v0.dot(v1);
  double d11 = v1.dot(v1);
  double d20 = v2.dot(v0);
  double d21 = v2.dot(v1);
  double denom = d00 * d11 - d01 * d01;
  v = (d11 * d20 - d01 * d21) / denom;
  w = (d00 * d21 - d01 * d20) / denom;
  u = 1.0 - v - w;
}

}

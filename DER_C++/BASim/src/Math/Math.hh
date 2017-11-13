/**
 * \file Math.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/03/2009
 */

#ifndef MATH_HH
#define MATH_HH

#include "BASim/src/Core/Util.hh"

namespace BASim {

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609;

/** Returns the square of the input value. */
template <typename T> inline T square(const T& x) { return x * x; }

/** Returns the cube of the input value. */
template <typename T> inline T cube(const T& x) { return x * x * x; }

/** Clamps scalar to the range [min,max]. */
template <typename T> inline T clamp( const T& scalar, const T& min, const T& max) 
{
  if (scalar < min) return min;
  if (scalar > max) return max;
  return scalar;
}

/** The Kronecker delta function. Returns 1 if \p i = \p j and 0
    otherwise. */
inline int KroneckerDelta(int i, int j)
{
  return (i == j ? 1 : 0);
}

/** Parallel transports vector \p u between two given tangent vectors
    \p t1 and \p t2. This corresponds to the rotation that aligns \p
    t1 to \p t2 while leaving unchanged all vectors orthogonal to the
    plane spanned by \p t1 and \p t2.

    \param[in] u The vector to be transported
    \param[in] t1 The first tangent (assumed to be normalized)
    \param[in] t2 The second tangent (assumed to be normalized)
    \return The parallel transported vector
*/
inline Vec3d parallelTransport(const Vec3d& u, Vec3d t1, Vec3d t2)
{
  Scalar d = t1.dot(t2);
  if (d == 1 || d == -1) return u;

  Vec3d b = t1.cross(t2);
  return d*u + 1/(1+d) * b.dot(u) * b + b.cross(u);
}

inline Vec3d parallel_transport(const Vec3d& u,
                                const Vec3d& t1,
                                const Vec3d& t2)
{
  Vec3d b = t1.cross(t2);
  if (b.norm() == 0) return u;
  b.normalize();
  b = (b - (b.dot(t1) * t1)).normalized();
  b = (b - (b.dot(t2) * t2)).normalized();
  Vec3d n1 = t1.cross(b);
  Vec3d n2 = t2.cross(b);
  return u.dot(t1) * t2 + u.dot(n1) * n2 + u.dot(b) * b;
}

/** Finds an orthonormal vector to the given input vector.

    \param[out] v The output vector that will be orthonormal to u
    \param[in] u The input vector
*/
inline void findOrthogonal(Vec3d& v, const Vec3d& u)
{
  assert(u.norm() != 0);

  v.setZero();
  int max = 0;
  for (int i = 0; i < u.size(); ++i) {
    if (u[i] == 0) {
      v[i] = 1;
      return;
    }
    if (fabs(u[i]) > fabs(u[max])) max = i;
  }

  int idx = (max + 1) % u.size();
  v[idx] = u[max];
  v[max] = -u[idx];
  v.normalize();

  assert(approxEq(u.dot(v), 0.0));
}

/** Computes the curvature binormal given two unit vectors.

    \param[out] curvatureBinormal The compute curvature binormal
    \param[in] t0 The first tangent
    \param[in] t1 The second tangent
*/
inline void computeCurvatureBinormal(Vec3d& curvatureBinormal,
                                     const Vec3d& t0,
                                     const Vec3d& t1)
{
  assert(approxEq(t0.norm(), 1.0));
  assert(approxEq(t1.norm(), 1.0));

  curvatureBinormal = 2.0 * t0.cross(t1) / (1.0 + t0.dot(t1));
}

/** Computes the curvature binormal given three vertices.

    \param[out] curvatureBinormal The computed curvature binormal
    \param[in] x0 The first vertex
    \param[in] x1 The second vertex
    \param[in] x2 The third vertex
*/
inline void computeCurvatureBinormal(Vec3d& curvatureBinormal,
                                     const Vec3d& x0,
                                     const Vec3d& x1,
                                     const Vec3d& x2)
{
  Vec3d t0 = (x1 - x0).normalized();
  Vec3d t1 = (x2 - x1).normalized();
  computeCurvatureBinormal(curvatureBinormal, t0, t1);
}

inline Scalar angle(const Vec3d& u, const Vec3d& v)
{
  return atan2(u.cross(v).norm(), u.dot(v));
}

/** Computes the signed angle from one vector to another given an
    orientation vector.

    \param[in] u The "from" vector.
    \param[in] v The "to" vector.
    \param[in] n The vector giving the orientation of the rotation plane.
    \return The signed angle between the two vectors
*/
inline Scalar signedAngle(const Vec3d& u, const Vec3d& v, const Vec3d& n)
{
  Vec3d w = u.cross(v);
  Scalar angle = atan2(w.norm(), u.dot(v));
  if (n.dot(w) < 0) return -angle;
  return angle;
}

/** Rotates a vector about the given axis (unit vector) by the given
    angle. The result is computed using Rodrigues' rotation formula
    (see <a
    href="http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula">here</a>).

    \param[in] v The vector to rotate. On return, it will contain the rotated vector.
    \param[in] z The (unit-length) axis about which to perform the rotation.
    \param[in] theta The angle of rotation.
*/
inline void rotateAxisAngle(Vec3d& v, Vec3d z, const Scalar& theta)
{
  assert(approxEq(z.norm(), 1.0));

  if (theta == 0) return;

  Scalar c = cos(theta);
  Scalar s = sin(theta);
  v = c * v + s * z.cross(v) + z.dot(v) * (1.0 - c) * z;
}

/** Returns the outer product of two vectors. In vector notation, this
    is \f$ab^T\f$. */
inline Mat3d outerProd(const Vec3d& a, const Vec3d& b)
{
  return a * b.transpose();
}

/** Returns the matrix representation of the cross product operator
    associated to a vector. Given a vector \f$a\f$it returns a matrix
    \f$A\f$ defined so that \f$Ab = a\times b\f$ for all vectors
    \f$b\f$. */
inline void crossMat(Mat3d& A, const Vec3d& a)
{
  A <<      0, -a(2),  a(1),
    a(2),     0, -a(0),
    -a(1),  a(0),     0;
}

/** \sa crossMat(Mat3d& A, const Vec3d& a) */
inline Mat3d crossMat(const Vec3d& a)
{
  Mat3d A;
  crossMat(A, a);
  return A;
}

/** Matrix representation of a quaternion. */
inline Mat4d quatToMat(const Quaternion& q)
{
  Mat4d Q;
  for (int i = 0; i < 4; ++i) Q(i,i) = q.w();
  Scalar a = q.x();
  Scalar b = q.y();
  Scalar c = q.z();
  Q(0,1) = -a; Q(1,0) =  a;
  Q(0,2) = -b; Q(2,0) =  b;
  Q(0,3) = -c; Q(3,0) =  c;
  Q(1,2) = -c; Q(2,1) =  c;
  Q(1,3) =  b; Q(3,1) = -b;
  Q(2,3) = -a; Q(3,2) =  a;
  return Q;
}

} // namespace BASim

#endif // MATH_HH

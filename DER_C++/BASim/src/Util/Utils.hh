#ifndef UTILS_HH
#define UTILS_HH

template <typename Derived> std::ostream&
operator<< (std::ostream& os, std::vector<Eigen::MatrixBase<Derived> >& v)
{
  os << "{";
  for (size_t i = 0; i < v.size(); ++i) {
    if (v[i].isVector()) os << v[i].format(EIGEN_VECTOR_IO);
    else os << v[i].format(EIGEN_MATRIX_IO);
    if (i < v.size() - 1) os << ",";
  }
  os << "}";
  return os;
}

template <class T> std::ostream& operator<<
(std::ostream& os, std::vector<T>& v)
{
  os << "{";
  for (size_t i = 0; i < v.size(); ++i) {
    os << v[i];
    if (i < v.size()-1) os << ",";
  }
  os << "}";
  return os;
}

namespace Utils {
// returns the rotation matrix representing a rotation about axis b
// by angle whose sine is s and cosine is c
inline void
matrixAxisSinCos(Mat3d& R, Vec3d b, const Scalar& s, const Scalar& c)
{
  if (s == 0 && c == 1) {
    R.setZero();
    R(0,0) = R(1,1) = R(2,2) = 1;
    return;
  }

  b.normalize();

  Scalar C = 1-c;
  Scalar xs = b[0]*s;   Scalar ys = b[1]*s;   Scalar zs = b[2]*s;
  Scalar xC = b[0]*C;   Scalar yC = b[1]*C;   Scalar zC = b[2]*C;
  Scalar xyC = b[0]*yC; Scalar yzC = b[1]*zC; Scalar zxC = b[2]*xC;

  R(0,0) = b[0]*xC+c;
  R(0,1) = xyC-zs;
  R(0,2) = zxC+ys;

  R(1,0) = xyC+zs;
  R(1,1) = b[1]*yC+c;
  R(1,2) = yzC-xs;

  R(2,0) = zxC-ys;
  R(2,1) = yzC+xs;
  R(2,2) = b[2]*zC+c;
}

// returns the rotation matrix representing a rotation of theta
// degrees about axis b
inline void matrixAxisAngle(Mat3d& M, const Vec3d& b, const Scalar& theta)
{
  matrixAxisSinCos(M, b, sin(theta), cos(theta));
}

inline void findOrthogonal(Mat3d& m, Vec3d v)
{
  assert(v.size() == m.cols()+1);
  assert(v.size() == m.rows());
  assert(v.norm() > 1.0e-6);

  v.normalize();

  int idx = 0;
  for (int i = 1; i < v.size(); i++) {
    if (fabs(v[i]) > fabs(v[idx])) idx = i;
  }

  int basis = 0;
  Vec3d w;
  for (int i = 0; i < v.size(); i++) {
    if (i == idx) continue;
    w.setZero();
    w[i] = v[idx];
    w[idx] = -v[i];
    w.normalize();
    m(i,basis) = w[i];
    m(idx,basis) = w[idx];
    basis++;
  }
}

/** Computes the holonomy of the closed loop formed by parallel
    transporting along \f$e^0 \to e^1\f$. */
inline Scalar computePsi(Vec3d e0, Vec3d e1, Vec3d e2, Vec3d e3)
{
  Vec3d u;
  Utils::findOrthogonal(e0, u);
  Vec3d u0(u);

  u = parallelTransport(u, e0.normalized(), e1.normalized());
  u = parallelTransport(u, e1.normalized(), e2.normalized());
  u = parallelTransport(u, e2.normalized(), e3.normalized());
  u = parallelTransport(u, e3.normalized(), e0.normalized());

  return signedAngle(u0, u, e0);
}

} // namespace Utils

#endif // UTILS_HH

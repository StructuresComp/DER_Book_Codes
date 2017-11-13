#include "BASim/src/Physics/ElasticRods/RodAnisoForce.hh"
#include "BASim/src/Math/Math.hh"

namespace BASim {

RodAnisoForce::RodAnisoForce(ElasticRod& rod)
  : RodAnisoBending(rod, "RodAnisoForce")
{}

Scalar RodAnisoForce::globalEnergy()
{
  Scalar E = 0.0;

  for (int i = 1; i < m_rod.nv() - 1; ++i) {

    Mat2d B = (getB(i - 1) + getB(i)) / 2.0;
    Scalar len = getRefVertexLength(i);

    //Scalar m = m_rod.twist(i);
    //E += _kt * m * m / (2.0 * len);

    for (int j = i - 1; j <= i; ++j) {

      const Vec2d& w = getOmega(i, j);
      const Vec2d& w0 = getOmegaBar(i, j);
      //const Mat2d& I = m_rod.I(j);

      E += 1.0 / (4.0 * len) * (w - w0).dot(B * (w - w0));
    }
  }

  return E;
}

void RodAnisoForce::globalForce(VecXd& F)
{
  _nk.resize(3*(m_rod.nv()-2));
  _np.resize(3*(m_rod.nv()-2));

  forPsi.resize(2*(m_rod.nv()-2));
  forPsi.setZero();

  computeNablaPsi();
  computeNablaKappa();

  std::vector<Vec3d> m1(m_rod.ne());
  std::vector<Vec3d> m2(m_rod.ne());
  for (int j = 0; j < m_rod.ne(); ++j) {
    m1[j] = m_rod.getMaterial1(j);
    m2[j] = m_rod.getMaterial2(j);
  }

  Scalar sum = 0;
  for (int i = m_rod.nv() - 2; i >= 1; --i) {
    Scalar len = getRefVertexLength(i);
    Vec3d forKappa(Vec3d::Zero());
    Mat2d B = (getB(i - 1) + getB(i)) / 2.0;

    for (int j = i - 1; j <= i; ++j) {
      const Vec2d& w = getOmega(i,j);
      const Vec2d& w0 = getOmegaBar(i, j);

      Vec2d temp = 1.0 / (2.0 * len) * (w - w0).transpose() * B;
      forKappa += Vec3d( temp(0) * m2[j](0) - temp(1) * m1[j](0) ,
                         temp(0) * m2[j](1) - temp(1) * m1[j](1) ,
                         temp(0) * m2[j](2) - temp(1) * m1[j](2) );

      forPsi[i-1+j] = -temp.dot(J * w);
    }

    Scalar forP = sum + forPsi[i - 1 + i];
    sum += forPsi[i - 1 + i] + forPsi[i - 1 + i - 1];

    for (int v = i-1; v <= i+1; ++v) {
      Mat3d& nk = evalNablaKappa(i, v);
      Vec3d& np = evalNablaPsi(i, v);
      Vec3d f = forP * np + nk.transpose() * forKappa;
      for (int cmp = 0; cmp < 3; ++cmp) {
        F(m_rod.vertIdx(v, cmp)) -= f(cmp);
      }
    }
  }
  /* THIS IS FOR QUASISTATICS -- COMMENTED OUT Sep 11 2009
  // nabla theta^j contributions
  if (m_rod.edgeFixed(m_rod.ne() - 1)) {
    int i = m_rod.nv() - 2;
    int j = i;
    Scalar len = getRefVertexLength(i);
    Scalar m = m_rod.twist(i);
    const Vec2d& w = getOmega(i,j);
    const Vec2d& w0 = getOmegaBar(i,j);
    const Mat2d& I = m_rod.I(j);
    Scalar s = 1.0 / (2.0 * len) * ((w - w0).transpose() * I).dot(J * w)
      - _kt * m / len;

    for (int i = 1; i < m_rod.nv()-1; ++i) {
      for (int v = i-1; v <= i+1; ++v) {
        Vec3d& np = evalNablaPsi(i, v);

        for (int cmp = 0; cmp < 3; ++cmp) F(3*v+cmp) -= s * np(cmp);
      }
    }
    }*/

  /* THIS IS FOR NON QUASISTATICS -- ADDED Sep 11 2009 */
  for (int i = 1; i < m_rod.nv() - 1; ++i) {
    Mat2d B = (getB(i - 1) + getB(i)) / 2.0;
    Scalar len = getRefVertexLength(i);

    for (int j = i - 1; j <= i; ++j) {
      const Vec2d& w = getOmega(i,j);
      const Vec2d& w0 = getOmegaBar(i,j);

      F(m_rod.edgeIdx(j)) -= 1.0 / (2.0 * len)
        * (w.transpose() * J).dot(B * (w - w0));
    }
  }
}

void RodAnisoForce::globalJacobian(int baseidx, Scalar scale, MatrixBase& Jacobian)
{
  assert(!"RodAnisoForce::globalJacobian not implemented");
}

void RodAnisoForce::computeNablaKappa()
{
  for (int i = 1; i < m_rod.nv()-1; ++i) {

    const Vec3d& kb = m_rod.getCurvatureBinormal(i);
    const Vec3d& e0 = m_rod.getEdge(i-1);
    Scalar len0 = m_rod.getEdgeLength(i-1);
    const Vec3d& e1 = m_rod.getEdge(i);
    Scalar len1 = m_rod.getEdgeLength(i);

    Mat3d& nk0 = _nk[3*(i-1)+0];
    crossMat(nk0, e1);
    nk0 *= 2;
    nk0 += outerProd(kb,e1);
    nk0 *= 1.0/(len0*len1 + e0.dot(e1));

    Mat3d& nk2 = _nk[3*(i-1)+2];
    crossMat(nk2, e0);
    nk2 *= 2;
    nk2 -= outerProd(kb,e0);
    nk2 *= 1.0/(len0*len1 + e0.dot(e1));

    Mat3d& nk1 = _nk[3*(i-1)+1];
    nk1 = nk0+nk2;
    nk1 *= -1;
  }
}

Mat3d& RodAnisoForce::evalNablaKappa(int i, int v)
{
  return _nk[2*(i-1)+v];
}

// evaluates grad (kb)_i with respect to vertex v
void RodAnisoForce::evalNablaKappa(Mat3d& nk, int i, int v)
{
  assert(i > 0);
  assert(i < ((int) m_rod.nv())-1);
  assert(fabs(i-v) <= 1);

  nk = _nk[2*(i-1)+v];
}

// finite difference approximation of grad (kb)_i with respect to vertex v
void RodAnisoForce::fdNablaKappa(Mat3d& nk, int i, int v)
{
  Scalar h = 1.0e-4;
  Vec3d e0 = m_rod.getEdge(i-1); Scalar len0 = m_rod.getEdgeLength(i-1);
  Vec3d e1 = m_rod.getEdge(i); Scalar len1 = m_rod.getEdgeLength(i);

  Scalar dir0 = 0, dir1 = 0;
  if (v == i-1) { dir0 = -1; dir1 = 0; }
  else if (v == i+1) { dir0 = 0; dir1 = 1; }
  else if (i == v) { dir0 = 1; dir1 = -1; }
  else assert(0);

  for (int k = 0; k < 3; ++k) {

    e0[k] += dir0*h; e1[k] += dir1*h;
    Vec3d kbp = 2.0/(len0*len1+e0.dot(e1)) * e0.cross(e1);

    e0[k] -= 2*dir0*h; e1[k] -= 2*dir1*h;
    Vec3d kbm = 2.0/(len0*len1+e0.dot(e1)) * e0.cross(e1);

    e0[k] += dir0*h; e1[k] += dir1*h;

    for (int cmp = 0; cmp < 3; ++cmp)
      nk(cmp,k) = (kbp[cmp]-kbm[cmp])/(2*h);
  }
}

void RodAnisoForce::computeNablaPsi()
{
  for (int i = 1; i < m_rod.nv()-1; ++i) {

    const Vec3d& kb = m_rod.getCurvatureBinormal(i);

    Scalar len0 = m_rod.getEdgeLength(i-1);
    Scalar len1 = m_rod.getEdgeLength(i);

    Vec3d np0 = 0.5/len0 * kb;
    Vec3d np2 = -0.5/len1 * kb;

    _np[3*(i-1)+0] = np0;
    _np[3*(i-1)+1] = -np0-np2;
    _np[3*(i-1)+2] = np2;

  }
}

Vec3d& RodAnisoForce::evalNablaPsi(int i, int v)
{
  return _np[2*(i-1) + v];
}

// evaluates grad psi_i with respect to vertex v
void RodAnisoForce::evalNablaPsi(Vec3d& np, int i, int v)
{
  assert(i > 0);
  assert(i < ((int)m_rod.nv())-1);
  assert(fabs(i-v) <= 1);

  np += _np[2*(i-1) + v];
}
/*
void RodAnisoForce::evalTorque(Vec3d& t, int j)
{
  Scalar mag2 = 0;

  const Mat2d& I = m_rod.I(j);

  for (int i = std::max(j, 1); i <= std::min(j + 1, (int) m_rod.nv() - 2); ++i) {
    Scalar m = m_rod.twist(i);

    Scalar len = getRefVertexLength(i);
    const Vec2d& w = getOmega(i,j);
    const Vec2d& w0 = getOmegaBar(i,j);

    //mag2 -= _kb/len * dot( w-w0, I*J*w );
    mag2 -= _kb/len * ((w-w0).transpose()*I).dot(J*w);
    mag2 += 2*_kt/len * m;
  }
  if (j == 0) mag2 *= -1;

  // direction of the torque is along edge j
  Vec3d dir = m_rod.getEdge(j);
  dir.normalize();
  dir *= mag2;
  t -= dir;
}
*/

void RodAnisoForce::updateProperties()
{
  RodAnisoBending::updateProperties();
}

} // namespace BASim

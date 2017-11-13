/**
 * \file TrackBall.cc
 */

#include "TrackBall.hh"

namespace {

inline void
rotate(BASim::Scalar m[4][4], const BASim::Vec3d& axis,
       const BASim::Scalar radians)
{
  BASim::Vec3d naxis( axis );
  naxis.normalize();

  BASim::Scalar x = naxis[0];
  BASim::Scalar y = naxis[1];
  BASim::Scalar z = naxis[2];

  BASim::Scalar c = cos( radians );
  BASim::Scalar s = sin( radians );
  BASim::Scalar t = 1 - c;

  m[0][0] = t*x*x + c;
  m[0][1] = t*x*y - z*s;
  m[0][2] = t*x*z + y*s;
  m[1][0] = t*x*y + z*s;
  m[1][1] = t*y*y + c;
  m[1][2] = t*y*z - x*s;
  m[2][0] = t*x*z - y*s;
  m[2][1] = t*y*z + x*s;
  m[2][2] = t*z*z + c;

  m[0][3] = 0;
  m[1][3] = 0;
  m[2][3] = 0;
  m[3][0] = 0;
  m[3][1] = 0;
  m[3][2] = 0;
  m[3][3] = 1;
}

} // namespace (global)

namespace BASim {

TrackBall::TrackBall(Camera* c)
  : m_camera(c)
  , m_rotating(false)
  , m_mode(GIMBAL)
{}

void TrackBall::start(const Vec2d& p)
{
  assert(p[0] >= -1 && p[0] <= 1 && p[1] >= -1 && p[1] <= 1);
  m_rotating = true;
  m_startPos = p;
  for (int i = 0; i < 4; ++i)
    m_rotation[i] = 0;
}

void TrackBall::update(const Vec2d& p)
{
  if (!m_rotating)
    return;

  Scalar m[4][4];
  if (m_startPos != p) {
    const Scalar coef(M_PI/2.0);

    const Scalar left_right_motion = p[0] - m_startPos[0];
    const Scalar up_down_motion = p[1] - m_startPos[1];

    // rotate about the 'up' vector
    Vec3d up;
    if (m_mode == GIMBAL) up = Vec3d(0,1,0);
    else                  up = m_camera->getUp();
    int sign = (up.dot(m_camera->getUp()) > 0) ? 1 : -1;
    rotate(m, up, coef * sign * left_right_motion);
    m_camera->orbit(m);

    // rotate about the horizontal vector
    Vec3d horizontal =
      m_camera->getUp().cross(m_camera->getEye() - m_camera->getViewCenter());
    rotate(m, horizontal, -coef * up_down_motion);
    m_camera->orbit(m);

    m_startPos = p;
  }
}

void TrackBall::stop()
{
  m_rotating = false;
}

void TrackBall::getRotation(Scalar r[4]) const
{
  for (int i = 0; i < 4; ++i) {
    r[i] = m_rotation[i];
  }
}

TrackBall::TrackBallMode TrackBall::getTrackBallMode() const
{
  return m_mode;
}

void TrackBall::setTrackBallMode(const TrackBall::TrackBallMode m)
{
  m_mode = m;
}

} // namespace BASim

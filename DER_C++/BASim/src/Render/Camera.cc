/**
 * \file Camera.cc
 */

#include "Camera.hh"
#include "OpenGLHeaders.hh"

namespace BASim {

Camera::Camera()
  : m_viewCenter(0,0,-1)
  , m_eye(0,0,0)
  , m_up(0,1,0)
  , m_dir(0,0,-1)
  , m_projMode(Camera::PERSPECTIVE)
//  , m_projMode(Camera::ORTHOGRAPHIC)
  , m_projParams(2)
  , m_zClipping(2)
{
  m_projParams[0] = 54.43;
  m_projParams[1] = 1;
  m_zClipping[0] = 10;
  m_zClipping[1] = 300;
}

Camera::Camera(const Camera& other)
  : m_viewCenter(other.m_viewCenter)
  , m_eye(other.m_eye)
  , m_up(other.m_up)
  , m_dir(other.m_dir)
  , m_projMode(other.m_projMode)
  , m_projParams(other.m_projParams)
  , m_zClipping(other.m_zClipping)
{}

Camera& Camera::operator=(const Camera& other)
{
  m_viewCenter = other.m_viewCenter;
  m_eye = other.m_eye;
  m_up = other.m_up;
  m_dir = other.m_dir;
  m_projMode = other.m_projMode;
  m_projParams = other.m_projParams;
  m_zClipping = other.m_zClipping;
  return *this;
}

void Camera::setOrthographic(const Scalar left, const Scalar right,
                             const Scalar top, const Scalar bottom)
{
  m_projMode = ORTHOGRAPHIC;
  m_projParams.resize(4);
  m_projParams[0] = left;
  m_projParams[1] = right;
  m_projParams[2] = top;
  m_projParams[3] = bottom;
}

void Camera::setPerspective(const Scalar fovy, const Scalar aspect)
{
  m_projMode = PERSPECTIVE;
  m_projParams.resize(2);
  m_projParams[0] = fovy;
  m_projParams[1] = aspect;
}

const Vec3d& Camera::getViewCenter() const
{
  return m_viewCenter;
}

const Vec3d& Camera::getEye() const
{
  return m_eye;
}

const Vec3d& Camera::getUp() const
{
  return m_up;
}

void Camera::setViewCenter(const Vec3d& c)
{
  m_viewCenter = c;
  m_dir = (m_viewCenter - m_eye).normalized();
}

void Camera::setEye(const Vec3d& e)
{
  m_eye = e;
}

void Camera::setUp(const Vec3d& u)
{
  m_up = u;
}

void Camera::getViewParams(Camera::ProjMode* m, std::vector<Scalar>* p) const
{
  assert(m && p);
  *m = m_projMode;
  p->resize(m_projParams.size());
  std::copy(m_projParams.begin(), m_projParams.end(), p->begin());
}

void Camera::setViewport(const int width, const int height)
{
  m_viewport.resize(2);
  m_viewport[0] = width;
  m_viewport[1] = height;
}

void Camera::getViewport(int* width, int* height) const
{
  assert(m_viewport.size() == 2);
  if (width) *width = m_viewport[0];
  if (height) *height = m_viewport[1];
}

void Camera::setZClipping(const Scalar near, const Scalar far)
{
  assert(m_zClipping.size() == 2);
  m_zClipping[0] = std::min(m_zClipping[0], near);
  m_zClipping[1] = std::max(m_zClipping[1], far);
}

void Camera::getZClipping(Scalar& near, Scalar& far) const
{
  assert(m_zClipping.size() == 2);
  near = m_zClipping[0];
  far = m_zClipping[1];
}

void Camera::setDefault3D(const Scalar sceneRadius)
{
  m_dir = Vec3d(0, 1, 0).normalized();
  m_eye = m_viewCenter - sceneRadius * 3.0 * m_dir;
  //m_eye = m_viewCenter + sceneRadius * 3.0f * Vec3d(1,1,1).normalized();
  //m_eye = m_viewCenter + sceneRadius * 3.0f * Vec3d(0,0,1).normalized();
  m_up = Vec3d(1,0,0);
}

void Camera::setDefault2D(const Scalar sceneRadius)
{
  m_eye = m_viewCenter + sceneRadius * 3.0f * Vec3d(0,0,1).normalized();
  m_up = Vec3d(0,1,0);
}
    
void Camera::setExperimentCamera()
{
  m_eye = Vec3d(-52.5645, 28.430269, 20.195164);
  m_up = Vec3d(0.33977593, 0.91739009, -0.20723836);
  m_dir = Vec3d(60.5645, -28.430269, -26.555327);
  m_viewCenter = m_eye + m_dir;
  m_dir.normalize();
}

void Camera::translateEye(const Vec3d& v)
{
  m_eye += v;
}

void Camera::translateCenter(const Vec3d& v)
{
  m_viewCenter += v;
}

// The lame-o Scalar pointer representation matches some legacy code we use.
void Camera::orbit(const Scalar m[4][4])
{
  {
    const Vec3d e = m_eye - m_viewCenter;
    m_eye[0] = e[0] * m[0][0] + e[1] * m[1][0] + e[2] * m[2][0];
    m_eye[1] = e[0] * m[0][1] + e[1] * m[1][1] + e[2] * m[2][1];
    m_eye[2] = e[0] * m[0][2] + e[1] * m[1][2] + e[2] * m[2][2];
    m_eye += m_viewCenter;
  }

  {
    const Vec3d u = m_up;
    m_up[0] = u[0] * m[0][0] + u[1] * m[1][0] + u[2] * m[2][0];
    m_up[1] = u[0] * m[0][1] + u[1] * m[1][1] + u[2] * m[2][1];
    m_up[2] = u[0] * m[0][2] + u[1] * m[1][2] + u[2] * m[2][2];
  }

  {
    const Vec3d d = m_dir;
    m_dir[0] = d[0] * m[0][0] + d[1] * m[1][0] + d[2] * m[2][0];
    m_dir[1] = d[0] * m[0][1] + d[1] * m[1][1] + d[2] * m[2][1];
    m_dir[2] = d[0] * m[0][2] + d[1] * m[1][2] + d[2] * m[2][2];
  }
}

void Camera::getSpanningSet(Vec3d* r, Vec3d* u, Vec3d* c) const
{
  //const Vec3d n = (m_viewCenter - m_eye).normalized();
  const Vec3d n = m_dir;
  const Vec3d w2 = (m_up - m_up.dot(n) * n).normalized();
  const Vec3d w1 = n.cross(w2);

  if (r) *r = w1;
  if (u) *u = w2;
  if (c) *c = n;
}

void Camera::applyProjection() const
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (m_projMode == ORTHOGRAPHIC) {
    assert(m_projParams.size() == 4);
    glOrtho(m_projParams[0], m_projParams[1],
            m_projParams[2], m_projParams[3],
            m_zClipping[0], m_zClipping[1]);

  } else if (m_projMode == PERSPECTIVE) {
    assert(m_projParams.size() == 2);
    gluPerspective(m_projParams[0],
                   m_projParams[1] * m_viewport[0] / (GLfloat)m_viewport[1],
                   m_zClipping[0], m_zClipping[1]);

  } else {
    assert(!"Unknown mode");
  }

  glMatrixMode(GL_MODELVIEW);
}

void Camera::applyViewport() const
{
  glViewport(0, 0, m_viewport[0], m_viewport[1]);
}

void Camera::applyCamera() const
{
  Vec3d center = m_eye + m_dir;
  gluLookAt(m_eye[0], m_eye[1], m_eye[2],                       // Eye position
            center[0], center[1], center[2],              // Look-at position
            m_up[0], m_up[1], m_up[2]);                         // Up vector
}

} // namespace BASim

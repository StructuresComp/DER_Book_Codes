/**
 * \file Zoomer.cc
 */

#include "Core/EigenIncludes.hh"
#include "Core/STLIncludes.hh"
#include "Core/Definitions.hh"
#include "OpenGLHeaders.hh"
#include "Camera.hh"
#include "Zoomer.hh"

namespace BASim {

Zoomer::Zoomer(Camera* c, const Scalar s)
  : m_camera(c)
  , m_translating(false)
  , m_scale(s)
{}

void Zoomer::setCamera(Camera* c)
{
  m_camera = c;
}

void Zoomer::setScale(const Scalar s)
{
  m_scale = s;
}

void Zoomer::start(const Vec2d& p)
{
  m_translating = true;
  m_startPos = p;
}

void Zoomer::update(const Vec2d& p)
{
  if (!m_translating) return;

  assert(m_camera);
  Vec3d in;
  m_camera->getSpanningSet(NULL, NULL, &in);

  const Vec3d translation = in * m_scale * (p[1] - m_startPos[1]);

  m_camera->translateEye(translation);

  m_startPos = p;
}

void Zoomer::stop()
{
  m_translating = false;
}

} // namespace BASim

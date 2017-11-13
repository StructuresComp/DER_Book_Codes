/**
 * \file RenderUtils.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/22/2009
 */

#ifndef RENDERUTILS_HH
#define RENDERUTILS_HH

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Render/OpenGLDecl.hh"
#include "BASim/src/Math/Math.hh"

namespace BASim {

inline void drawArrow(const Vec3d& start, Vec3d dir, Scalar scale = 1.0)
{
  dir.normalize();
  Vec3d end = start + scale * dir;
  Scalar len = (end - start).norm();
  Scalar rad = 0.03 * len;

  Vec3d x;
  findOrthogonal(x, end - start);
  Vec3d y = (end - start).cross(x).normalized();

  glBegin(GL_QUAD_STRIP);
  for (unsigned int i = 0; i < 33; ++i) {
    Scalar angle = 2.0 * M_PI * i / 32.0;
    Vec3d norm = cos(angle) * x - sin(angle) * y;
    Vec3d point = rad * norm;

    OpenGL::normal(norm);
    OpenGL::vertex((start + point).eval());

    point = start + 0.8 * (end - start) + point;
    OpenGL::normal(norm);
    OpenGL::vertex(point);
  }
  glEnd();

  Vec3d center = start + 0.8 * (end - start);
  Vec3d orientation = (end - start).normalized();
  Vec3d axis = Vec3d(0,0,1).cross(orientation);
  Scalar angle = 180.0/M_PI*atan2(axis.norm(), Vec3d(0,0,1).dot(orientation));
  if (angle != 0) axis.normalize();
  glPushMatrix();
  glTranslated(center(0), center(1), center(2));
  glRotated(angle, axis(0), axis(1), axis(2));
  glutSolidCone(2.2 * rad, 0.2 * len, 32, 2);
  glPopMatrix();
}

} // namespace BASim

#endif // RENDERUTILS_HH

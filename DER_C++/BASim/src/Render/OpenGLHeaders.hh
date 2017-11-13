/**
 * \file OpenGLHeaders.hh
 */

#ifndef OPENGLHEADERS_HH
#define OPENGLHEADERS_HH

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#  include <GLUT/glut.h>
#else
#  include <GL/gl.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

#endif // OPENGLHEADERS_HH

/**
 * \file RenderBase.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/11/2009
 */

#include "BASim/src/Render/RenderBase.hh"

namespace BASim {

void RenderBase::render()
{
  if (enableNormals())  setNormalArray();
  if (enableColors())   setColorArray();
  if (enableVertices()) drawVertices();
  if (enableNormals())  unsetNormalArray();
  if (enableColors())   unsetColorArray();
}

void RenderBase::drawVertices()
{
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_COLOR_MATERIAL);
  glColor3ub(255, 255, 255);

  // activate and specify pointer to vertex array
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_FLOAT, 0, &(vertices[0]));

  // draw the triangles
  glDrawElements(GL_TRIANGLES, vertices.size() / 3, GL_UNSIGNED_INT, &(indices[0]));

  // deactivate vertex array
  glDisableClientState(GL_VERTEX_ARRAY);

  glDisable(GL_COLOR_MATERIAL);
}

void RenderBase::setNormalArray()
{
  glEnableClientState(GL_NORMAL_ARRAY);
  glNormalPointer(GL_FLOAT, 0, &(normals[0]));
}

void RenderBase::unsetNormalArray()
{
  glDisableClientState(GL_NORMAL_ARRAY);
}

void RenderBase::setColorArray()
{
  glEnableClientState(GL_COLOR_ARRAY);
  glColorPointer(3, GL_UNSIGNED_BYTE, 0, &(colors[0]));
}

void RenderBase::unsetColorArray()
{
  glDisableClientState(GL_COLOR_ARRAY);
}

} // namespace BASim

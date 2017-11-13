/**
 * \file Scene.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 08/30/2009
 */

#include "Scene.hh"

namespace BASim {

void Scene::render()
{
  for (size_t i = 0; i < m_renderers.size(); ++i) {
    m_renderers[i]->render();
  }
}

void Scene::addRenderer(RenderBase& renderer)
{
  m_renderers.push_back(&renderer);
}

} // namespace BASim

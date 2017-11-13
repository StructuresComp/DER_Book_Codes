/**
 * \file Scene.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/30/2009
 */

#ifndef SCENE_HH
#define SCENE_HH

#include "BASim/src/Render/RenderBase.hh"

namespace BASim {

/** Represents a scene to be rendered. */
class Scene
{
public:

  Scene() {}

  ~Scene() {}

  void render();

  void addRenderer(RenderBase& renderer);

protected:

  std::vector<RenderBase*> m_renderers;
};

} // namespace BASim

#endif // SCENE_HH

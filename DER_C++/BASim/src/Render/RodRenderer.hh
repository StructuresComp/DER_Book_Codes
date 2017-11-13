/**
 * \file RodRenderer.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/30/2009
 */

#ifndef RODRENDERER_HH
#define RODRENDERER_HH

#include "BASim/src/Render/RenderBase.hh"
#include "BASim/src/Render/Color.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Physics/ElasticRods/RodTube.hh"

namespace BASim {

/** Class that implements OpenGL rendering for rods. */
class RodRenderer : public RenderBase
{
public:

  enum DrawMode { SIMPLE, SMOOTH, PHOTOREALISTIC, NONE };

  RodRenderer(ElasticRod& rod);

  void render();

  DrawMode getMode() const { return m_mode; }
  void setMode(DrawMode mode) { m_mode = mode; }

  bool& drawMaterial() { return m_drawMaterial; }
  bool& drawReference() { return m_drawReference; }
  bool& scaleToRadius() { return m_scaleToRadius; }
  bool& drawArrows() { return m_drawArrows; }

  virtual Vec3d calculateObjectCenter();
  virtual Scalar calculateObjectBoundingRadius(const Vec3d& center);

  static RodRenderer * lastInstance() { return s_last_instance; }
    
protected:

  void drawSimpleRod();
  void drawSmoothRod();
  void drawPRRod();
    
  void drawMaterialFrame();
  void drawReferenceFrame();

  ElasticRod& m_rod;
  RodTube m_tube;

  DrawMode m_mode;

  bool m_drawMaterial;
  bool m_drawReference;
  bool m_scaleToRadius;
  bool m_drawArrows;

  std::vector<Color> m_palette;
    
  static RodRenderer * s_last_instance;
    
};

} // namespace BASim

#endif // RODRENDERER_HH

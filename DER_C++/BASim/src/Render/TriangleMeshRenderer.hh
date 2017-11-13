/**
 * \file TriangleMeshRenderer.hh
 *
 * \author smith@cs.columbia.edu
 * \date 03/15/2010
 */

#ifndef TRIANGLEMESHRENDERER_HH
#define TRIANGLEMESHRENDERER_HH

#include "BASim/src/Core/TriangleMesh.hh"
#include "BASim/src/Render/RenderBase.hh"

namespace BASim {

  /** Class that implements OpenGL rendering for triangle meshes. */
  class TriangleMeshRenderer : public RenderBase
  {
  public:
  
    enum DrawMode { DBG, FLAT, NONE };

    TriangleMeshRenderer( const TriangleMesh& mesh );
    
    void render();
    
    DrawMode getMode() const { return m_mode; }
    void setMode(DrawMode mode) { m_mode = mode; }
    
    virtual Vec3d calculateObjectCenter();
    virtual Scalar calculateObjectBoundingRadius(const Vec3d& center);
    
    static TriangleMeshRenderer * lastInstance() { return s_last_instance; }
    const TriangleMesh & mesh() const { return m_mesh; }
      
  protected:
    const TriangleMesh& m_mesh;
    DrawMode m_mode;
      
    static TriangleMeshRenderer * s_last_instance;
      
  };
  
} // namespace BASim

#endif // TRIANGLEMESHRENDERER_HH

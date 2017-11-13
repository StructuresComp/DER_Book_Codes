/**
 * \file TriangleMeshRenderer.cc
 *
 * \author smith@cs.columbia.edu
 * \date 03/15/2010
 */

#include "TriangleMeshRenderer.hh"
#include "BASim/src/Render/Color.hh"
#include "BASim/src/Render/OpenGLDecl.hh"
#include "RodRenderer.hh"

namespace BASim 
{

TriangleMeshRenderer * TriangleMeshRenderer::s_last_instance = NULL;
    
TriangleMeshRenderer::TriangleMeshRenderer( const TriangleMesh& mesh )
: m_mesh(mesh)
, m_mode(FLAT)
{
    s_last_instance = this;
}
  
void TriangleMeshRenderer::render()
{
  // Render all edges
  glLineWidth(2);
  glBegin(GL_LINES);
  OpenGL::color(Color(0,0,0));
  for( TriangleMesh::edge_iter eit = m_mesh.edges_begin(); eit != m_mesh.edges_end(); ++eit )
  {
    OpenGL::vertex(m_mesh.getVertex(m_mesh.fromVertex(*eit)));
    OpenGL::vertex(m_mesh.getVertex(m_mesh.toVertex(*eit)));
  }
  glEnd();
  
  // Render all faces
  glBegin(GL_TRIANGLES);
  OpenGL::color(Color(130,130,130));
  for( TriangleMesh::face_iter fit = m_mesh.faces_begin(); fit != m_mesh.faces_end(); ++fit )
  {
    for( TriangleMesh::FaceVertexIter fvit = m_mesh.fv_iter(*fit); fvit; ++fvit )
    {
//      OpenGL::vertex(m_mesh.getVertex(*fvit));
    }      
  }
  glEnd();
  
  // Render all vertices
  glPointSize(5);
  glBegin(GL_POINTS);
  OpenGL::color(Color(0,0,0));
  for( TriangleMesh::vertex_iter vit = m_mesh.vertices_begin(); vit != m_mesh.vertices_end(); ++vit ) 
    OpenGL::vertex(m_mesh.getVertex(*vit));
  glEnd();
  
  
  {
		TriangleMesh::vertex_iter vit = m_mesh.vertices_begin();
		
		Vec3d x0 = m_mesh.getVertex(*vit); ++vit;
		Vec3d x1 = m_mesh.getVertex(*vit); ++vit;
		Vec3d x2 = m_mesh.getVertex(*vit); ++vit;
    
		Scalar t = x0[0];
		Scalar b = x2[0];
		Scalar r = -12;
		Scalar l = 12;

      if (!RodRenderer::lastInstance() || RodRenderer::lastInstance()->getMode() != RodRenderer::PHOTOREALISTIC)
      {
		glLineWidth(1);
//		glColor3ub(0, 0, 0);
  OpenGL::color(Color(0,0,0));
  
		while (t < b) {
			glBegin(GL_LINES);
glColor3ub(0, 0, 0);
  OpenGL::color(Color(0,0,0));
				glVertex3f(t, 0.5, r);
				glVertex3f(t, 0.5, r+5.0);
			glEnd();

			glBegin(GL_LINES);
glColor3ub(0, 0, 0);
  OpenGL::color(Color(0,0,0));
				glVertex3f(t, 0.5, l-5.0);
				glVertex3f(t, 0.5, l);
			glEnd();

			t += 5.0;
		}
      }
	}
	
}

Vec3d TriangleMeshRenderer::calculateObjectCenter()
{
  Vec3d center(0.0,0.0,0.0);

  return center;
  
  for( TriangleMesh::vertex_iter vit = m_mesh.vertices_begin(); vit != m_mesh.vertices_end(); ++vit ) 
  {
    center += m_mesh.getVertex(*vit);
  }

  if( m_mesh.nv() != 0 ) center /= ((double)m_mesh.nv());
  
  return center;
}

double TriangleMeshRenderer::calculateObjectBoundingRadius( const Vec3d& center )
{
  Scalar radius = 0.0;
  
  return radius;
    
  for( TriangleMesh::vertex_iter vit = m_mesh.vertices_begin(); vit != m_mesh.vertices_end(); ++vit )
  {
    radius = std::max(radius, (m_mesh.getVertex(*vit) - center).norm());
  }
  
  return radius;
}

} // namespace BASim

/**
 * \file TriangleMesh.hh
 *
 * \author smith@cs.columbia.edu
 * \date 03/15/2010
 */

#ifndef TRIANGLEMESH_HH
#define TRIANGLEMESH_HH

#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"

namespace BASim {

/**
 * A triangle mesh.
 */
class TriangleMesh : public TopologicalObject
{
public:
  
  TriangleMesh()
  {
    add_property(m_vertex_positions, "vertex_positions", Vec3d(0,0,0));
    add_property(m_vertex_velocities, "vertex_velocities", Vec3d(0,0,0));
    fixed = false;
  }
  
  /**
   * Access and modify the mesh's vertices.
   */
  Vec3d& getVertex( const TopologicalObject::vertex_handle& vh )
  {
    return property(m_vertex_positions)[vh];
  }

  Vec3d& getVelocity( const TopologicalObject::vertex_handle& vh )
  {
    return property(m_vertex_velocities)[vh];
  }

  void setVelocity( const Vec3d& new_vel )
  {
		for(vertex_iter vIt = vertices_begin(); vIt != vertices_end(); ++vIt) {
			getVelocity(*vIt) = new_vel;
		}
  }
  
  /**
   * Const version.
   */
  const Vec3d& getVertex( const TopologicalObject::vertex_handle& vh ) const
  {
    return property(m_vertex_positions)[vh];
  }

  const Vec3d& getVelocity( const TopologicalObject::vertex_handle& vh ) const
  {
    return property(m_vertex_velocities)[vh];
  }  
  
  bool fixed;

private:
  VPropHandle<Vec3d> m_vertex_positions;
  VPropHandle<Vec3d> m_vertex_velocities;
  
};

} // namespace BASim

#endif // TRIANGLEMESH_HH

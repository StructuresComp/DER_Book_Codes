/**
 * \file TopologicalObject.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/12/2009
 * \author smith@cs.columbia.edu
 * \date 03/14/2010
 */

#ifndef TOPOLOGICALOBJECT_HH
#define TOPOLOGICALOBJECT_HH

#include "BASim/src/Core/ObjectBase.hh"
#include "BASim/src/Core/TopologicalObject/TopObjHandles.hh"
#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"
#include "BASim/src/Core/TopologicalObject/Topology.hh"

namespace BASim {

/** An object that represents a collection of vertices, edges and faces
   with associated connectivity information. Warning: Face support is
   still largely untested and subject to revision. */
class TopologicalObject : public ObjectBase
{
public:

  typedef VertexHandle<TopologicalObject>               vertex_handle;
  typedef VertexIterator<TopologicalObject>             vertex_iter;
  typedef const VertexIterator<TopologicalObject>       const_vertex_iter;
  typedef VertexTopology<TopologicalObject>             vertex_topology;

  typedef EdgeHandle<TopologicalObject>                 edge_handle;
  typedef EdgeIterator<TopologicalObject>               edge_iter;
  typedef const EdgeIterator<TopologicalObject>         const_edge_iter;
  typedef EdgeTopology<TopologicalObject>               edge_topology;

  typedef FaceHandle<TopologicalObject>                 face_handle;
  typedef FaceIterator<TopologicalObject>               face_iter;
  typedef const FaceIterator<TopologicalObject>         const_face_iter;
  typedef FaceTopology<TopologicalObject>               face_topology;

  typedef VertexEdgeIterator<TopologicalObject>         VertexEdgeIter;
  typedef const VertexEdgeIterator<TopologicalObject>   ConstVertexEdgeIter;

  typedef EdgeVertexIterator<TopologicalObject>         EdgeVertexIter;
  typedef const EdgeVertexIterator<TopologicalObject>   ConstEdgeVertexIter;

  typedef VertexVertexIterator<TopologicalObject>       VertexVertexIter;
  typedef const VertexVertexIterator<TopologicalObject> ConstVertexVertexIter;

  typedef FaceVertexIterator<TopologicalObject>         FaceVertexIter;
  typedef const FaceVertexIterator<TopologicalObject>   ConstFaceVertexIter;

  TopologicalObject();

  virtual ~TopologicalObject() {}

  /** Number of vertices */
  int nv() const { return property(m_nv); }
  /** Number of edges */
  int ne() const { return property(m_ne); }
  /** Number of faces */
  int nf() const { return property(m_nf); }

  /** Add a vertex to the object */
  vertex_handle addVertex();
  /** Add a directed edge that connects the two given vertices */
  edge_handle addEdge(const vertex_handle& v0, const vertex_handle& v1);
  /** Add a triangular face to the object */
  face_handle addFace(const vertex_handle& v0, const vertex_handle& v1, const vertex_handle& v2);

  /** Deletes the given vertex and associated edges, faces (NOT IMPLEMENTED) */
  void deleteVertex(const vertex_handle& vertex);
  /** Deletes the given edge and associated faces (NOT IMPLEMENTED) */
  void deleteEdge(const edge_handle& edge);
  /** Deletes the given face (NOT IMPLEMENTED) */
  void deleteFace(const face_handle& face);

  /** Returns a handle to the vertex at the head of the edge */
  vertex_handle fromVertex(const edge_handle& eh) const;
  /** Returns a handle to the vertex at the tail of the edge */
  vertex_handle toVertex(const edge_handle& eh) const;

  /** \name Iterators */

  //@{

  /** Returns an iterator to the first vertex */
  vertex_iter vertices_begin();
  /** Const version */
  const_vertex_iter vertices_begin() const;

  /** Returns an iterator to the past-the-end vertex */
  vertex_iter vertices_end();
  /** Const version */
  const_vertex_iter vertices_end() const;

  /** Returns an iterator to the first edge */
  edge_iter edges_begin();
  /** Const version */
  const_edge_iter edges_begin() const;

  /** Returns an iterator to the past-the-end edge */
  edge_iter edges_end();
  /** Const version */
  const_edge_iter edges_end() const;

  /** Returns an iterator to the first face */
  face_iter faces_begin();
  /** Const version */
  const_face_iter faces_begin() const;

  /** Returns an iterator to the past-the-end face */
  face_iter faces_end();
  /** Const version */
  const_face_iter faces_end() const;

  //@}

  /** \name Circulators */

  //@{

  /** Circulator over edges adjacent to a vertex */
  VertexEdgeIter ve_iter(const vertex_handle& vh);
  /** Const version */
  ConstVertexEdgeIter ve_iter(const vertex_handle& vh) const;

  /** Circulator over vertices adjacent to an edge */
  EdgeVertexIter ev_iter(const edge_handle& eh);
  /** Const version */
  ConstEdgeVertexIter ev_iter(const edge_handle& eh) const;

  /** Circulator over vertices in the one-ring of the given vertex */
  VertexVertexIter vv_iter(const vertex_handle& vh);
  /** Const version */
  ConstVertexVertexIter vv_iter(const vertex_handle& vh) const;

  /** Circulator over vertices adjacent to a face */
  FaceVertexIter fv_iter(const face_handle& fh);
  /** Const version */
  ConstFaceVertexIter fv_iter(const face_handle& fh) const;

  //@}

  /** \name Topology */

  //@{

  /** Returns topology associated to a vertex */
  vertex_topology& getVertexTopology(const vertex_handle& vh);
  /** Const version */
  const vertex_topology& getVertexTopology(const vertex_handle& vh) const;

  /** Returns topology associated to an edge */
  edge_topology& getEdgeTopology(const edge_handle& eh);
  /** Const version */
  const edge_topology& getEdgeTopology(const edge_handle& eh) const;

  /** Returns topology associated to a face */
  face_topology& getFaceTopology(const face_handle& fh);
  /** Const version */
  const face_topology& getFaceTopology(const face_handle& eh) const;

  //@}

  /** \name Vertex properties */
  //@{
  BA_CREATE_PROPERTY(VPropHandle, m_vertexProps, nv());
  //@}

  /** \name Edge properties */
  //@{
  BA_CREATE_PROPERTY(EPropHandle, m_edgeProps, ne());
  //@}

  /** \name Face properties */
  //@{
  BA_CREATE_PROPERTY(FPropHandle, m_faceProps, nf());
  //@}

  BA_INHERIT_BASE(ObjectBase);

protected:

  PropertyContainer m_vertexProps; ///< Vertex property container
  PropertyContainer m_edgeProps;   ///< Edge property container
  PropertyContainer m_faceProps;   ///< Face property container

  ObjPropHandle<int> m_nv; ///< number of vertices
  ObjPropHandle<int> m_ne; ///< number of edges
  ObjPropHandle<int> m_nf; ///< number of faces

  VPropHandle<vertex_topology> m_vertTop; ///< Topology of the vertices
  EPropHandle<edge_topology> m_edgeTop;   ///< Topology of the edges
  FPropHandle<face_topology> m_faceTop;   ///< Topology of the faces
};

#include "TopologicalObject.inl"

} // namespace BASim

#endif // TOPOLOGICALOBJECT_HH

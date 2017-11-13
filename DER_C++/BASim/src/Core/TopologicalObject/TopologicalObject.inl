/**
 * \file TopologicalObject.inl
 *
 * \author miklos@cs.columbia.edu
 * \date 09/12/2009
 * \author smith@cs.columbia.edu
 * \date 03/14/2010
 */

inline TopologicalObject::TopologicalObject()
{
  add_property(m_nv, "number_of_vertices", 0);
  add_property(m_ne, "number_of_edges", 0);
  add_property(m_nf, "number_of_faces", 0);

  add_property(m_vertTop, "vertex_topology");
  add_property(m_edgeTop, "edge_topology");
  add_property(m_faceTop, "face_topology");
}

inline TopologicalObject::vertex_handle TopologicalObject::addVertex()
{
  int idx = nv();
  property(m_nv) += 1;
  m_vertexProps.resize(nv());

  return vertex_handle(idx);
}

inline TopologicalObject::edge_handle
TopologicalObject::addEdge(const vertex_handle& v0, const vertex_handle& v1)
{
  int idx = ne();
  property(m_ne) += 1;
  m_edgeProps.resize(ne());
  edge_handle ehnd(idx);

  edge_topology& etop = getEdgeTopology(ehnd);
  etop.setFromVertex(v0);
  etop.setToVertex(v1);

  vertex_topology& vtop0 = getVertexTopology(v0);
  vtop0.addEdge(ehnd);

  vertex_topology& vtop1 = getVertexTopology(v1);
  vtop1.addEdge(ehnd);

  return ehnd;
}

inline TopologicalObject::face_handle
TopologicalObject::addFace(const vertex_handle& v0, const vertex_handle& v1, const vertex_handle& v2)
{
  // TODO: Add checks that vertex handles are valid
  
  int idx = nf();
  property(m_nf) += 1;
  m_faceProps.resize(nf());
  face_handle fhnd(idx);
  
  face_topology& ftop = getFaceTopology(fhnd);
  assert( ftop.size() == 3 );
  ftop.setVertex(0,v0);
  ftop.setVertex(1,v1);
  ftop.setVertex(2,v2);

  // TODO: add face information to VertexTopology about faces that it belongs to

  return fhnd;
}

inline void TopologicalObject::deleteVertex(const vertex_handle& vertex)
{
//  bool delete_vertex_is_implemnted = false;
//  assert( delete_vertex_is_implemnted );

//	std::cout << vertex.idx() << " delete vertex in topObj.inl\n";
	
  vertex_topology& vtop = getVertexTopology(vertex);
  
  while(vtop.size() > 0) {
  	deleteEdge(vtop[0]);
  }

  property(m_nv) -= 1;
  m_vertexProps.delete_element(vertex);
  
	for(edge_iter eIt = edges_begin(); eIt != edges_end(); ++eIt) {
		getEdgeTopology(*eIt).adjustVertexIdxAfterDeleteVertex(vertex);
	}
  
}

inline void TopologicalObject::deleteEdge(const edge_handle& edge)
{
//  bool delete_edge_is_implemnted = false;
//  assert( delete_edge_is_implemnted );

//	std::cout << edge.idx() << " delete edge in topObj.inl\n";

  edge_topology& etop = getEdgeTopology(edge);
  
  getVertexTopology(etop.getFromVertex()).deleteEdge(edge);
  getVertexTopology(etop.getToVertex()).deleteEdge(edge);
  
  property(m_ne) -= 1;
  m_edgeProps.delete_element(edge);

	for(vertex_iter vIt = vertices_begin(); vIt != vertices_end(); ++vIt) {
		getVertexTopology(*vIt).adjustEdgeIdxAfterDeleteEdge(edge);
	}

}

inline void TopologicalObject::deleteFace(const face_handle& face)
{
  bool delete_face_is_implemnted = false;
  if( !delete_face_is_implemnted ) std::cout << "DELETE FACE NOT SUPPORTED" << std::endl;
  assert( delete_face_is_implemnted );
}

inline TopologicalObject::vertex_iter 
TopologicalObject::vertices_begin()
{
  return vertex_iter(this, vertex_handle(0));
}

inline TopologicalObject::const_vertex_iter
TopologicalObject::vertices_begin() const
{
  return vertex_iter(this, vertex_handle(0));
}

inline TopologicalObject::vertex_iter 
TopologicalObject::vertices_end()
{
  return vertex_iter(this, vertex_handle(nv()));
}

inline TopologicalObject::const_vertex_iter
TopologicalObject::vertices_end() const
{
  return vertex_iter(this, vertex_handle(nv()));
}

inline TopologicalObject::edge_iter 
TopologicalObject::edges_begin()
{
  return edge_iter(this, edge_handle(0));
}

inline TopologicalObject::const_edge_iter 
TopologicalObject::edges_begin() const
{
  return edge_iter(this, edge_handle(0));
}

inline TopologicalObject::face_iter 
TopologicalObject::faces_begin()
{
  return face_iter(this, face_handle(0));
}

inline TopologicalObject::const_face_iter 
TopologicalObject::faces_begin() const
{
  return face_iter(this, face_handle(0));
}

inline TopologicalObject::edge_iter 
TopologicalObject::edges_end()
{
  return edge_iter(this, edge_handle(ne()));
}

inline TopologicalObject::const_edge_iter 
TopologicalObject::edges_end() const
{
  return edge_iter(this, edge_handle(ne()));
}

inline TopologicalObject::face_iter 
TopologicalObject::faces_end()
{
  return face_iter(this, face_handle(nf()));
}

inline TopologicalObject::const_face_iter 
TopologicalObject::faces_end() const
{
  return face_iter(this, face_handle(nf()));
}

inline TopologicalObject::vertex_topology&
TopologicalObject::getVertexTopology(const vertex_handle& vh)
{
  return property(m_vertTop)[vh];
}

inline const TopologicalObject::vertex_topology&
TopologicalObject::getVertexTopology(const vertex_handle& vh) const
{
  return property(m_vertTop)[vh];
}

inline TopologicalObject::edge_topology&
TopologicalObject::getEdgeTopology(const edge_handle& eh)
{
  return property(m_edgeTop)[eh];
}

inline const TopologicalObject::edge_topology&
TopologicalObject::getEdgeTopology(const edge_handle& eh) const
{
  return property(m_edgeTop)[eh];
}

inline TopologicalObject::face_topology&
TopologicalObject::getFaceTopology(const face_handle& fh)
{
  return property(m_faceTop)[fh];
}

inline const TopologicalObject::face_topology&
TopologicalObject::getFaceTopology(const face_handle& fh) const
{
  return property(m_faceTop)[fh];
}

inline TopologicalObject::VertexEdgeIter
TopologicalObject::ve_iter(const vertex_handle& vh)
{
  return VertexEdgeIter(this, vh);
}

inline TopologicalObject::ConstVertexEdgeIter
TopologicalObject::ve_iter(const vertex_handle& vh) const
{
  return VertexEdgeIter(this, vh);
}

inline TopologicalObject::EdgeVertexIter
TopologicalObject::ev_iter(const edge_handle& eh)
{
  return EdgeVertexIter(this, eh);
}

inline TopologicalObject::ConstEdgeVertexIter
TopologicalObject::ev_iter(const edge_handle& eh) const
{
  return EdgeVertexIter(this, eh);
}

inline TopologicalObject::VertexVertexIter
TopologicalObject::vv_iter(const vertex_handle& vh)
{
  return VertexVertexIter(this, vh);
}

inline TopologicalObject::ConstVertexVertexIter
TopologicalObject::vv_iter(const vertex_handle& vh) const
{
  return VertexVertexIter(this, vh);
}

inline TopologicalObject::vertex_handle
TopologicalObject::fromVertex(const edge_handle& eh) const
{
  return getEdgeTopology(eh).getFromVertex();
}

inline TopologicalObject::vertex_handle
TopologicalObject::toVertex(const edge_handle& eh) const
{
  return getEdgeTopology(eh).getToVertex();
}

inline TopologicalObject::FaceVertexIter
TopologicalObject::fv_iter(const face_handle& fh)
{
  return FaceVertexIter(this, fh);
}

inline TopologicalObject::ConstFaceVertexIter
TopologicalObject::fv_iter(const face_handle& fh) const
{
  return FaceVertexIter(this, fh);
}

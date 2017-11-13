/**
 * \file VertexStencil.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/13/2009
 */

#ifndef VERTEXSTENCIL_HH
#define VERTEXSTENCIL_HH

#include "BASim/src/Physics/ElasticRods/Stencil.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"

namespace BASim {

class VertexStencil;

template <> struct TypeInfo<VertexStencil>
{
  typedef ElasticRod::vertex_handle handle_type;
  typedef ElasticRod::vertex_iter   iterator;
};

/** Class for iterating over the vertices of a rod and getting indices
    associated to degrees of freedom. */
class VertexStencil : public StencilT<VertexStencil>
{
public:

  VertexStencil(ElasticRod& obj)
    : StencilT<VertexStencil>(obj)
  {}

  virtual ~VertexStencil() {}

  virtual VertexStencil& operator++ ()
  {
    ++m_iter;
    return *this;
  }

  virtual handle_ref handle()
  {
    return *m_iter;
  }

  const_handle_ref handle() const
  {
    return *m_iter;
  }

  virtual iterator begin() const
  {
    return ++(m_obj->vertices_begin());
  }

  virtual iterator end() const
  {
    return --(m_obj->vertices_end());
  }

  virtual VertexStencil& operator= (const iterator& it)
  {
    m_iter = it;
    return *this;
  }

  virtual bool operator== (const iterator& it) const
  {
    return m_iter == it;
  }

  virtual bool operator!= (const iterator& it) const
  {
    return !(operator==(it));
  }

  virtual void indices(IndexArray& indices)
  {
    indices.resize(11);
    handle_ref h = handle();
    ElasticRod& rod = *smart_cast<ElasticRod*>(m_obj);
    ElasticRod::VertexVertexIter vvit = rod.vv_iter(h);
    handle_type h_prev = *vvit;
    ++vvit;
    handle_type h_next = *vvit;

    for (int i = 0; i < 3; ++i) {
      indices(i) = rod.vertIdx(h_prev, i);
      indices(4 + i) = rod.vertIdx(h, i);
      indices(8 + i) = rod.vertIdx(h_next, i);
    }

    ElasticRod::VertexEdgeIter veit = rod.ve_iter(h);

    indices(3) = rod.edgeIdx(*veit);
    ++veit;
    indices(7) = rod.edgeIdx(*veit);
  }

  ElasticRod::edge_handle inEdge() const
  {
    return getRod().getVertexTopology(handle())[0];
  }

  ElasticRod::edge_handle outEdge() const
  {
    return getRod().getVertexTopology(handle())[1];
  }

  ElasticRod::vertex_handle prevVert() const
  {
    return getRod().getEdgeTopology(inEdge())[0];
  }

  ElasticRod::vertex_handle nextVert() const
  {
    return getRod().getEdgeTopology(outEdge())[1];
  }

protected:

  inline const ElasticRod& getRod() const
  { return *smart_cast<const ElasticRod*>(m_obj); }
  inline ElasticRod& getRod() { return *smart_cast<ElasticRod*>(m_obj); }

  iterator m_iter;
};

} // namespace BASim

#endif // VERTEXSTENCIL_HH

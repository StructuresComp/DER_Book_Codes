/**
 * \file EdgeStencil.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/13/2009
 */

#ifndef EDGESTENCIL_HH
#define EDGESTENCIL_HH

#include "BASim/src/Physics/ElasticRods/Stencil.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"

namespace BASim {

class EdgeStencil;

template <> struct TypeInfo<EdgeStencil>
{
  typedef ElasticRod::edge_handle handle_type;
  typedef ElasticRod::edge_iter   iterator;
};

/** Class for iterating over the edges of a rod and getting indices
    associated to degrees of freedom. */
class EdgeStencil : public StencilT<EdgeStencil>
{
public:

  EdgeStencil(ElasticRod& obj)
    : StencilT<EdgeStencil>(obj)
  {}

  virtual ~EdgeStencil() {}

  virtual EdgeStencil& operator++ ()
  {
    ++m_iter;
    return *this;
  }

  virtual handle_ref handle()
  {
    return *m_iter;
  }

  virtual const_handle_ref handle() const
  {
    return *m_iter;
  }

  virtual iterator begin() const
  {
    return m_obj->edges_begin();
  }

  virtual iterator end() const
  {
    return m_obj->edges_end();
  }

  virtual EdgeStencil& operator= (const iterator& it)
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
    indices.resize(6);
    handle_ref h = handle();
    ElasticRod& rod = *smart_cast<ElasticRod*>(m_obj);
    for (int i = 0; i < 3; ++i) {
      indices(i) = rod.vertIdx(rod.fromVertex(h), i);
      indices(3 + i) = rod.vertIdx(rod.toVertex(h), i);
    }
  }

protected:

  iterator m_iter;
};

} // namespace BASim

#endif // EDGESTENCIL_HH

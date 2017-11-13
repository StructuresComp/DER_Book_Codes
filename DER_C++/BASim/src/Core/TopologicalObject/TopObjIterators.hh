/**
 * \file TopObjIterators.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/13/2009
 * \author smith@cs.columbia.edu
 * \date 03/14/2010
 */

#ifndef TOPOBJITERATORS_HH
#define TOPOBJITERATORS_HH

#include "BASim/src/Core/Definitions.hh"

namespace BASim {

/** Iterator over vertices */
template <class T>
class VertexIterator
{
public:

  typedef typename T::vertex_handle   value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  VertexIterator() : m_obj(NULL), m_hnd(-1) {}

  VertexIterator(const VertexIterator& vit)
    : m_obj(vit.m_obj)
    , m_hnd(vit.m_hnd)
  {}

  VertexIterator(const T* t, value_type val)
    : m_obj(t)
    , m_hnd(val)
  {}

  VertexIterator& operator= (const VertexIterator& vit)
  {
    m_obj = vit.m_obj;
    m_hnd = vit.m_hnd;
    return *this;
  }

  bool operator== (const VertexIterator& vit) const
  {
    return ((m_obj == vit.m_obj) && (m_hnd == vit.m_hnd));
  }

  bool operator!= (const VertexIterator& vit) const
  {
    return !(operator==(vit));
  }

  VertexIterator& operator++ ()
  {
    m_hnd.__increment();
    return *this;
  }

  VertexIterator& operator-- ()
  {
    m_hnd.__decrement();
    return *this;
  }

  reference operator* ()
  {
    return m_hnd;
  }

  const_reference operator* () const
  {
    return m_hnd;
  }

  pointer operator-> ()
  {
    return &m_hnd;
  }

  const_pointer operator-> () const
  {
    return &m_hnd;
  }

protected:

  const T* m_obj;
  value_type m_hnd;
};

/** Iterator over edges */
template <class T>
class EdgeIterator
{
public:

  typedef typename T::edge_handle   value_type;
  typedef value_type&               reference;
  typedef const value_type&         const_reference;
  typedef value_type*               pointer;
  typedef const value_type*         const_pointer;

  EdgeIterator() : m_obj(NULL), m_hnd(-1) {}

  EdgeIterator(const EdgeIterator& eit)
    : m_obj(eit.m_obj)
    , m_hnd(eit.m_hnd)
  {}

  EdgeIterator(const T* t, value_type val)
    : m_obj(t)
    , m_hnd(val)
  {}

  EdgeIterator& operator= (const EdgeIterator& eit)
  {
    m_obj = eit.m_obj;
    m_hnd = eit.m_hnd;
    return *this;
  }

  bool operator== (const EdgeIterator& eit) const
  {
    return ((m_obj == eit.m_obj) && (m_hnd == eit.m_hnd));
  }

  bool operator!= (const EdgeIterator& eit) const
  {
    return !(operator==(eit));
  }

  EdgeIterator& operator++ ()
  {
    m_hnd.__increment();
    return *this;
  }

  EdgeIterator& operator-- ()
  {
    m_hnd.__decrement();
    return *this;
  }

  reference operator* ()
  {
    return m_hnd;
  }

  const_reference operator* () const
  {
    return m_hnd;
  }

  pointer operator-> ()
  {
    return &m_hnd;
  }

  const_pointer operator-> () const
  {
    return &m_hnd;
  }

protected:

  const T* m_obj;
  value_type m_hnd;
};

/** Iterator over faces */
template <class T>
class FaceIterator
{
public:

  typedef typename T::face_handle   value_type;
  typedef value_type&               reference;
  typedef const value_type&         const_reference;
  typedef value_type*               pointer;
  typedef const value_type*         const_pointer;

  FaceIterator() : m_obj(NULL), m_hnd(-1) {}

  FaceIterator(const FaceIterator& fit)
    : m_obj(fit.m_obj)
    , m_hnd(fit.m_hnd)
  {}

  FaceIterator(const T* t, value_type val)
    : m_obj(t)
    , m_hnd(val)
  {}

  FaceIterator& operator= (const FaceIterator& fit)
  {
    m_obj = fit.m_obj;
    m_hnd = fit.m_hnd;
    return *this;
  }

  bool operator== (const FaceIterator& fit) const
  {
    return ((m_obj == fit.m_obj) && (m_hnd == fit.m_hnd));
  }

  bool operator!= (const FaceIterator& fit) const
  {
    return !(operator==(fit));
  }

  FaceIterator& operator++ ()
  {
    m_hnd.__increment();
    return *this;
  }

  FaceIterator& operator-- ()
  {
    m_hnd.__decrement();
    return *this;
  }

  reference operator* ()
  {
    return m_hnd;
  }

  const_reference operator* () const
  {
    return m_hnd;
  }

  pointer operator-> ()
  {
    return &m_hnd;
  }

  const_pointer operator-> () const
  {
    return &m_hnd;
  }

protected:

  const T* m_obj;
  value_type m_hnd;
};

/** Circulator over edges adjacent to a vertex */
template <class T>
class VertexEdgeIterator
{
public:

  typedef typename T::vertex_topology topology;
  typedef typename T::edge_handle     edge_handle;
  typedef typename T::vertex_handle   vertex_handle;
  typedef edge_handle                 value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  VertexEdgeIterator() : m_obj(NULL), m_hnd(-1), m_top(NULL), m_idx(0) {}

  VertexEdgeIterator(const VertexEdgeIterator& veit)
    : m_obj(veit.m_obj)
    , m_hnd(veit.m_hnd)
    , m_top(veit.m_top)
    , m_idx(0)
  {}

  VertexEdgeIterator(const T* t, const vertex_handle& vh)
    : m_obj(t)
    , m_hnd(vh)
    , m_top(&t->getVertexTopology(vh))
    , m_idx(0)
  {}

  VertexEdgeIterator& operator= (const VertexEdgeIterator& veit)
  {
    m_obj = veit.m_obj;
    m_hnd = veit.m_hnd;
    m_top = veit.m_top;
    m_idx = veit.m_idx;
    return *this;
  }

  bool operator== (const VertexEdgeIterator& veit) const
  {
    assert(m_top);
    return ((m_obj == veit.m_obj)
            && ((*m_top)[m_idx] == (*veit.m_top)[veit.m_idx]));
  }

  bool operator!= (const VertexEdgeIterator& veit) const
  {
    return !(operator==(veit));
  }

  VertexEdgeIterator& operator++ ()
  {
    ++m_idx;
    return *this;
  }

  VertexEdgeIterator& operator-- ()
  {
    --m_idx;
    return *this;
  }

  reference operator* ()
  {
    assert(m_top);
    return (*const_cast<topology*>(m_top))[m_idx];
  }

  const_reference operator* () const
  {
    assert(m_top);
    return (*m_top)[m_idx];
  }

  pointer operator-> ()
  {
    assert(m_top);
    return &(*m_top)[m_idx];
  }

  const_pointer operator-> () const
  {
    assert(m_top);
    return &(*m_top)[m_idx];
  }

  operator bool() const
  {
    assert(m_top);
    return m_idx < (*m_top).size();
  }

protected:

  const T* m_obj;        ///< Pointer to the topological object
  vertex_handle m_hnd;   ///< Handle to the vertex whose topology is being considered
  const topology* m_top; ///< Pointer to connectivity information for vertex
  size_t m_idx;          ///< Current index into topological information
};

/** Circulator over vertices adjacent to an edge */
template <class T>
class EdgeVertexIterator
{
public:

  typedef typename T::edge_topology   topology;
  typedef typename T::edge_handle     edge_handle;
  typedef typename T::vertex_handle   vertex_handle;
  typedef vertex_handle               value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  EdgeVertexIterator() : m_obj(NULL), m_hnd(-1), m_top(NULL), m_idx(0) {}

  EdgeVertexIterator(const EdgeVertexIterator& evit)
    : m_obj(evit.m_obj)
    , m_hnd(evit.m_hnd)
    , m_top(evit.m_top)
    , m_idx(0)
  {}

  EdgeVertexIterator(const T* t, const edge_handle& eh)
    : m_obj(t)
    , m_hnd(eh)
    , m_top(&t->getEdgeTopology(eh))
    , m_idx(0)
  {}

  EdgeVertexIterator& operator= (const EdgeVertexIterator& evit)
  {
    m_obj = evit.m_obj;
    m_hnd = evit.m_hnd;
    m_top = evit.m_top;
    m_idx = evit.m_idx;
    return *this;
  }

  bool operator== (const EdgeVertexIterator& evit) const
  {
    assert(m_top);
    assert(evit.m_top);
    return ((m_obj == evit.m_obj)
            && ((*m_top)[m_idx] == (*evit.m_top)[evit.m_idx]));
  }

  bool operator!= (const EdgeVertexIterator& evit) const
  {
    return !(operator==(evit));
  }

  EdgeVertexIterator& operator++ ()
  {
    ++m_idx;
    return *this;
  }

  EdgeVertexIterator& operator-- ()
  {
    --m_idx;
    return *this;
  }

  reference operator* ()
  {
    assert(m_top);
    return (*const_cast<topology*>(m_top))[m_idx];
  }

  const_reference operator* () const
  {
    assert(m_top);
    return (*m_top)[m_idx];
  }

  pointer operator-> ()
  {
    assert(m_top);
    return &(*m_top)[m_idx];
  }

  const_pointer operator-> () const
  {
    assert(m_top);
    return &(*m_top)[m_idx];
  }

  operator bool() const
  {
    assert(m_top);
    return m_idx < (*m_top).size();
  }

protected:

  const T* m_obj;        ///< Pointer to topological object
  edge_handle m_hnd;     ///< Handle to edge whose connectivity is being considered
  const topology* m_top; ///< Pointer to connectivity information for edge
  size_t m_idx;          ///< Current index into topology information
};

/** Circulator over vertices in the one ring of a vertex */
template <class T>
class VertexVertexIterator
{
public:

  typedef typename T::vertex_topology topology;
  typedef typename T::edge_handle     edge_handle;
  typedef typename T::vertex_handle   vertex_handle;
  typedef vertex_handle               value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  VertexVertexIterator() : m_obj(NULL), m_hnd(-1), m_top(NULL), m_idx(0) {}

  VertexVertexIterator(const VertexVertexIterator& vvit)
    : m_obj(vvit.m_obj)
    , m_hnd(vvit.m_hnd)
    , m_top(vvit.m_top)
    , m_idx(0)
  {}

  VertexVertexIterator(const T* t, const vertex_handle& vh)
    : m_obj(t)
    , m_hnd(vh)
    , m_top(&t->getVertexTopology(vh))
    , m_idx(0)
  {}

  VertexVertexIterator& operator= (const VertexVertexIterator& vvit)
  {
    m_obj = vvit.m_obj;
    m_hnd = vvit.m_hnd;
    m_top = vvit.m_top;
    m_idx = vvit.m_idx;
    return *this;
  }

  bool operator== (const VertexVertexIterator& vvit) const
  {
    assert(m_top);
    assert(vvit.mtop);
    return ((m_obj == vvit.m_obj)
            && ((*m_top)[m_idx] == (*vvit.m_top)[vvit.m_idx]));
  }

  bool operator!= (const VertexVertexIterator& vvit) const
  {
    return !(operator==(vvit));
  }

  VertexVertexIterator& operator++ ()
  {
    ++m_idx;
    return *this;
  }

  VertexVertexIterator& operator-- ()
  {
    --m_idx;
    return *this;
  }

  reference operator* ()
  {
    assert(m_top);
    edge_handle& eh = (*const_cast<topology*>(m_top))[m_idx];
    m_target = m_obj->fromVertex(eh);
    if (m_target == m_hnd) m_target = m_obj->toVertex(eh);
    return m_target;
  }

  const_reference operator* () const
  {
    assert(m_top);
    const edge_handle& eh = (*m_top)[m_idx];
    m_target = m_obj->fromVertex(eh);
    if (m_target == m_hnd) m_target = m_obj->toVertex(eh);
    return m_target;
  }

  pointer operator-> ()
  {
    return &(this->operator*());
  }

  const_pointer operator-> () const
  {
    return &(this->operator*());
  }

  operator bool() const
  {
    assert(m_top);
    return m_idx < (*m_top).size();
  }

protected:

  const T* m_obj;         ///< Pointer to topological object
  vertex_handle m_hnd;    ///< Handle to vertex whose one-ring is being considered
  vertex_handle m_target; ///< Handle to vertex that dereference oprator returns
  const topology* m_top;  ///< Pointer to connectivity information for edge
  size_t m_idx;           ///< Current index into topology information
};

/** Circulator over vertices adjacent to a face */
template <class T>
class FaceVertexIterator
{
public:

  typedef typename T::face_topology   topology;
  typedef typename T::face_handle     face_handle;
  typedef typename T::vertex_handle   vertex_handle;
  typedef vertex_handle               value_type;
  typedef value_type&                 reference;
  typedef const value_type&           const_reference;
  typedef value_type*                 pointer;
  typedef const value_type*           const_pointer;

  FaceVertexIterator() : m_obj(NULL), m_hnd(-1), m_top(NULL), m_idx(0) {}

  FaceVertexIterator(const FaceVertexIterator& vfit)
    : m_obj(vfit.m_obj)
    , m_hnd(vfit.m_hnd)
    , m_top(vfit.m_top)
    , m_idx(0)
  {}

  FaceVertexIterator(const T* t, const face_handle& fh)
    : m_obj(t)
    , m_hnd(fh)
    , m_top(&t->getFaceTopology(fh))
    , m_idx(0)
  {}

  FaceVertexIterator& operator= (const FaceVertexIterator& vfit)
  {
    m_obj = vfit.m_obj;
    m_hnd = vfit.m_hnd;
    m_top = vfit.m_top;
    m_idx = vfit.m_idx;
    return *this;
  }

  bool operator== (const FaceVertexIterator& vfit) const
  {
    assert(m_top);
    assert(vfit.m_top);
    return ((m_obj == vfit.m_obj)
            && ((*m_top)[m_idx] == (*vfit.m_top)[vfit.m_idx]));
  }

  bool operator!= (const FaceVertexIterator& vfit) const
  {
    return !(operator==(vfit));
  }

  FaceVertexIterator& operator++ ()
  {
    ++m_idx;
    return *this;
  }

  FaceVertexIterator& operator-- ()
  {
    --m_idx;
    return *this;
  }

  reference operator* ()
  {
    assert(m_top);
    return (*const_cast<topology*>(m_top))[m_idx];
  }

  const_reference operator* () const
  {
    assert(m_top);
    return (*m_top)[m_idx];
  }

  pointer operator-> ()
  {
    assert(m_top);
    return &(*m_top)[m_idx];
  }

  const_pointer operator-> () const
  {
    assert(m_top);
    return &(*m_top)[m_idx];
  }

  operator bool() const
  {
    assert(m_top);
    return m_idx < (*m_top).size();
  }

protected:

  const T* m_obj;        ///< Pointer to topological object
  face_handle m_hnd;     ///< Handle to face whose connectivity is being considered
  const topology* m_top; ///< Pointer to connectivity information for face
  size_t m_idx;          ///< Current index into topology information
};

} // namespace BASim

#endif // TOPOBJITERATORS_HH

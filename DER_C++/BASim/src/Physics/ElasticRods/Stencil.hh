/**
 * \file Stencil.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/13/2009
 */

#ifndef STENCIL_HH
#define STENCIL_HH

#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"

namespace BASim {

/** Stencil interface. A stencil is essentially an iterator over an
    object with convenience methods for accessing the properties of
    the object associated to the stencil */
class Stencil
{
public:

  Stencil(TopologicalObject& obj) : m_obj(&obj) {}

  virtual ~Stencil() {}

protected:

  TopologicalObject* m_obj;
};

template <typename T> struct TypeInfo {};

template <class Derived>
class StencilT : public Stencil
{
public:

  typedef typename TypeInfo<Derived>::handle_type      handle_type;
  typedef handle_type&                                 handle_ref;
  typedef const handle_type&                           const_handle_ref;
  typedef handle_type*                                 handle_ptr;
  typedef const handle_type*                           const_handle_ptr;
  typedef typename TypeInfo<Derived>::iterator         iterator;
  typedef const typename TypeInfo<Derived>::iterator   const_iterator;

  StencilT(TopologicalObject& obj) : Stencil(obj) {}

  virtual Derived& operator++ () = 0;
  virtual handle_ref handle() = 0;
  virtual const_handle_ref handle() const = 0;
  virtual iterator begin() const = 0;
  virtual iterator end() const = 0;

  virtual Derived& operator= (const iterator& it) = 0;
  virtual bool operator== (const iterator& it) const = 0;
  virtual bool operator!= (const iterator& it) const = 0;
  virtual void indices(IndexArray& indices) = 0;
};

} // namespace BASim

#endif // STENCIL_HH

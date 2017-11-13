/**
 * \file TopObjHandles.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/13/2009
 * \author smith@cs.columbia.edu
 * \date 03/14/2010
 */

#ifndef TOPOBJHANDLES_HH
#define TOPOBJHANDLES_HH

#include "BASim/src/Core/Handle.hh"
#include "BASim/src/Core/Property.hh"

namespace BASim {

/** Handle for referring to a vertex */
template <class T>
class VertexHandle : public HandleBase
{
public:

  explicit VertexHandle(int idx = -1) : HandleBase(idx) {}
};

/** Handle for referring to an edge */
template <class T>
class EdgeHandle : public HandleBase
{
public:

  explicit EdgeHandle(int idx = -1) : HandleBase(idx) {}
};

/** Handle for referring to a face */
template <class T>
class FaceHandle : public HandleBase
{
public:
  explicit FaceHandle(int idx = -1) : HandleBase(idx) {}
};

/** Handle for referring to vertex properties. */
template <typename T>
class VPropHandle : public PropertyHandleBase<T>
{
public:

  explicit VPropHandle(int idx = -1) : PropertyHandleBase<T>(idx) {}
  explicit VPropHandle(const PropertyHandleBase<T>& b)
    : PropertyHandleBase<T>(b)
  {}
};

/** Handle for referring to edge properties. */
template <typename T>
class EPropHandle : public PropertyHandleBase<T>
{
public:

  explicit EPropHandle(int idx = -1) : PropertyHandleBase<T>(idx) {}
  explicit EPropHandle(const PropertyHandleBase<T>& b)
    : PropertyHandleBase<T>(b)
  {}
};

/** Handle for referring to face properties. */
template <typename T>
class FPropHandle : public PropertyHandleBase<T>
{
public:

  explicit FPropHandle(int idx = -1) : PropertyHandleBase<T>(idx) {}
  explicit FPropHandle(const PropertyHandleBase<T>& b)
    : PropertyHandleBase<T>(b)
  {}
};

} // namespace BASim

#endif // TOPOBJHANDLES_HH

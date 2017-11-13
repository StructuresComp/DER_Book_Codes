/**
 * \file ObjectHandles.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/13/2009
 */

#ifndef OBJECTHANDLES_HH
#define OBJECTHANDLES_HH

#include "BASim/src/Core/Handle.hh"
#include "BASim/src/Core/Property.hh"

namespace BASim {

/** Handle for referring to objects. */
class ObjectHandle : public HandleBase
{
public:

  explicit ObjectHandle(int idx = -1) : HandleBase(idx) {}
};

/** Handle for referring to object properties. */
template <typename T>
class ObjPropHandle : public PropertyHandleBase<T>
{
public:

  explicit ObjPropHandle(int idx = -1) : PropertyHandleBase<T>(idx) {}
  explicit ObjPropHandle(const PropertyHandleBase<T>& b)
    : PropertyHandleBase<T>(b)
  {}
};

} // namespace BASim

#endif // OBJECTHANDLES_HH

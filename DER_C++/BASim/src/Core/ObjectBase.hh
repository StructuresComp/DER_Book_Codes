/**
 * \file ObjectBase.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef OBJECTBASE_HH
#define OBJECTBASE_HH

#include "BASim/src/Core/Macros.hh"
#include "BASim/src/Core/Property.hh"
#include "BASim/src/Core/ObjectHandles.hh"

namespace BASim {

/** Base class for all objects. Objects derived from this class can
    have different properties added to them that will be associated to
    the object (rather than, for example, associated to vertices). */
class ObjectBase
{
public:

  ObjectBase() {}
  virtual ~ObjectBase() {}

  /** \name Property access

      Provides access to object properties and the
      PropertyContainer. */

  //@{

  BA_CREATE_SINGULAR_PROPERTY(ObjPropHandle, m_objectProps);

  //@}

protected:

  PropertyContainer m_objectProps; ///< properties associated to the object
};

} // namespace BASim

#endif // OBJECTBASE_HH

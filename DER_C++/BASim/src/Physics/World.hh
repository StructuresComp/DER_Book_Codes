/**
 * \file World.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/09/2009
 */

#ifndef WORLD_HH
#define WORLD_HH

#include "BASim/src/Core/ObjectBase.hh"
#include "BASim/src/Core/ObjectControllerBase.hh"

namespace BASim {

/** Class that collects all of the objects in the world. */
class World : public ObjectBase
{
public:

  typedef std::vector<ObjectBase*> Objects;
  typedef std::vector<ObjectControllerBase*> Controllers;

  World() {}
  ~World() {}

  void initialize(int argc, char** argv)
  {
#ifdef HAVE_PETSC
    PetscUtils::initializePetsc(&argc, &argv);
#endif // HAVE_PETSC
  }

  void finalize()
  {
#ifdef HAVE_PETSC
    PetscUtils::finalizePetsc();
#endif // HAVE_PETSC
  }

  ObjectHandle addObject(ObjectBase* object)
  {
    assert( object != NULL );

    int idx = m_objects.size();
    m_objects.push_back(object);

    return ObjectHandle(idx);
  }

  ObjectBase& getObject(const ObjectHandle& oh)
  {
    assert(oh.isValid());
    assert(oh.idx() >= 0);
    assert((size_t) oh.idx() < m_objects.size());

    return *m_objects[oh.idx()];
  }

  const Objects& getObjects() const
  {
    return m_objects;
  }

  Objects& getObjects()
  {
    return m_objects;
  }

  ObjectControllerHandle addController(ObjectControllerBase* controller)
  {
    assert( controller != NULL );
    
    int idx = m_controllers.size();
    m_controllers.push_back(controller);

    return ObjectControllerHandle(idx);
  }

  ObjectControllerBase& getController(const ObjectControllerHandle& och)
  {
    assert(och.isValid());
    assert(och.idx() >= 0);
    assert(och.idx() < (int) m_controllers.size());

    return *m_controllers[och.idx()];
  }

  Controllers& getControllers()
  {
    return m_controllers;
  }

  void execute()
  {
    Controllers::iterator it;
    for (it = m_controllers.begin(); it != m_controllers.end(); ++it) {
      (*it)->execute();
    }
  }
  
  void removeObject(ObjectBase* object)
  {
    Objects::iterator it;
    for (it = m_objects.begin(); it != m_objects.end(); ++it) {
      if ((*it) == object) {
        m_objects.erase(it);
        break;
      }
    }  
  }
  
  void clear()
  {
    m_objects.clear();
    m_controllers.clear();
  
  }

protected:

  Objects m_objects;
  Controllers m_controllers;
};

} // namespace BASim

#endif // WORLD_HH

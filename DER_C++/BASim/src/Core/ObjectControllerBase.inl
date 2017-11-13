/**
 * \file ObjectControllerBase.inl
 *
 * \author miklos@cs.columbia.edu
 * \date 09/05/2009
 */

inline ObjectHandle ObjectControllerBase::addObject(ObjectBase& object)
{
  int idx = m_objects.size();
  m_objects.push_back(&object);

  return ObjectHandle(idx);
}

inline ObjectBase&
ObjectControllerBase::getObject(const ObjectHandle& handle)
{
  assert(handle.isValid());
  assert(handle.idx() < (int) m_objects.size());

  return *m_objects[handle.idx()];
}

inline const ObjectBase&
ObjectControllerBase::getObject(const ObjectHandle& handle) const
{
  assert(handle.isValid());
  assert(handle.idx() < (int) m_objects.size());

  return *m_objects[handle.idx()];
}

inline std::vector<ObjectBase*>& ObjectControllerBase::getObjects()
{
  return m_objects;
}

inline const std::vector<ObjectBase*>&
ObjectControllerBase::getObjects() const
{
  return m_objects;
}

/**
 * \file PhysObject.inl
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

inline PhysObject::PhysObject()
{
  add_property(m_ndof, "number of degrees of freedom");
  add_property(m_map, "index map");
}

inline int PhysObject::ndof() const
{
  return property(m_ndof);
}

inline const Scalar& PhysObject::getDof(int i) const
{
  const DofHandle& handle = property(m_map).getDof(i);

  if (handle.getType() == DofHandle::VERTEX_DOF) {
    const vertex_handle& vh
      = static_cast<const vertex_handle&>(handle.getHandle());
    return getVertexDof(vh, handle.getNum());

  } else {
    const edge_handle& eh
      = static_cast<const edge_handle&>(handle.getHandle());
    return getEdgeDof(eh, handle.getNum());

  }
}

inline void PhysObject::setDof(int i, const Scalar& dof)
{
  const DofHandle& handle = property(m_map).getDof(i);

  if (handle.getType() == DofHandle::VERTEX_DOF) {
    const vertex_handle& vh
      = static_cast<const vertex_handle&>(handle.getHandle());
    setVertexDof(vh, handle.getNum(), dof);

  } else {
    const edge_handle& eh
      = static_cast<const edge_handle&>(handle.getHandle());
    setEdgeDof(eh, handle.getNum(), dof);

  }
}

inline const Scalar& PhysObject::getVel(int i) const
{
  const DofHandle& handle = property(m_map).getDof(i);

  if (handle.getType() == DofHandle::VERTEX_DOF) {
    const vertex_handle& vh
      = static_cast<const vertex_handle&>(handle.getHandle());
    return getVertexVel(vh, handle.getNum());

  } else {
    const edge_handle& eh
      = static_cast<const edge_handle&>(handle.getHandle());
    return getEdgeVel(eh, handle.getNum());

  }
}

inline void PhysObject::setVel(int i, const Scalar& dof)
{
  const DofHandle& handle = property(m_map).getDof(i);

  if (handle.getType() == DofHandle::VERTEX_DOF) {
    const vertex_handle& vh
      = static_cast<const vertex_handle&>(handle.getHandle());
    setVertexVel(vh, handle.getNum(), dof);

  } else {
    const edge_handle& eh
      = static_cast<const edge_handle&>(handle.getHandle());
    setEdgeVel(eh, handle.getNum(), dof);

  }
}

inline const Scalar& PhysObject::getMass(int i) const
{
  const DofHandle& handle = property(m_map).getDof(i);

  if (handle.getType() == DofHandle::VERTEX_DOF) {
    const vertex_handle& vh
      = static_cast<const vertex_handle&>(handle.getHandle());
    return getVertexMass(vh, handle.getNum());

  } else {
    const edge_handle& eh
      = static_cast<const edge_handle&>(handle.getHandle());
    return getEdgeMass(eh, handle.getNum());

  }
}

inline void PhysObject::setMass(int i, const Scalar& dof)
{
  const DofHandle& handle = property(m_map).getDof(i);

  if (handle.getType() == DofHandle::VERTEX_DOF) {
    const vertex_handle& vh
      = static_cast<const vertex_handle&>(handle.getHandle());
    setVertexMass(vh, handle.getNum(), dof);

  } else {
    const edge_handle& eh
      = static_cast<const edge_handle&>(handle.getHandle());
    setEdgeMass(eh, handle.getNum(), dof);

  }
}

/**
 * \file PhysObject.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef PHYSOBJECT_HH
#define PHYSOBJECT_HH

#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"
#include "BASim/src/Physics/DegreeOfFreedom.hh"

namespace BASim {

class MatrixBase;

/** Base class for physics objects. Provides an array-like interface
    to the degrees of freedom of the object. Derived classes and/or
    other classes that use a class derived from this one are
    responsible for setting up the indexing of the degrees of freedom
    (how to map indices to the degrees of freedom and vice-versa).
*/
class PhysObject : public TopologicalObject
{
public:

  PhysObject();

  virtual ~PhysObject() {}

  /** The number of degrees of freedom this object has. */
  int ndof() const;

  /** This function must be called after all initialization to the
      object has been made. */
  virtual void setup() {}

  /** \name Accessors

      Provide access to degrees of freedom, velocities, and masses. */

  //@{

  virtual const Scalar& getDof(int i) const;
  virtual void setDof(int i, const Scalar& dof);

  virtual const Scalar& getVel(int i) const;
  virtual void setVel(int i, const Scalar& vel);

  virtual const Scalar& getMass(int i) const;
  virtual void setMass(int i, const Scalar& mass);

  //@}

  /** Getting and setting vertex, edge degrees of freedom, velocities,
      and masses */

  //@{

  virtual const Scalar&
  getVertexDof(const vertex_handle& vh, int num) const = 0;

  virtual void
  setVertexDof(const vertex_handle& vh, int num, const Scalar& dof) = 0;

  virtual const Scalar&
  getEdgeDof(const edge_handle& eh, int num) const = 0;

  virtual void
  setEdgeDof(const edge_handle& eh, int num, const Scalar& dof) = 0;

  virtual const Scalar&
  getVertexVel(const vertex_handle& vh, int num) const = 0;

  virtual void
  setVertexVel(const vertex_handle& vh, int num, const Scalar& vel) = 0;

  virtual const Scalar&
  getEdgeVel(const edge_handle& eh, int num) const = 0;

  virtual void
  setEdgeVel(const edge_handle& eh, int num, const Scalar& vel) = 0;

  virtual const Scalar&
  getVertexMass(const vertex_handle& vh, int num) const = 0;

  virtual void
  setVertexMass(const vertex_handle& vh, int num, const Scalar& mass) = 0;

  virtual const Scalar&
  getEdgeMass(const edge_handle& eh, int num) const = 0;

  virtual void
  setEdgeMass(const edge_handle& eh, int num, const Scalar& mass) = 0;

  //@}

  /** \name Internal force and Jacobian computation */

  //@{

  virtual void computeForces(VecXd& force) {}
  virtual void computeJacobian(MatrixBase& J) {}

  //@}

protected:

  ObjPropHandle<int> m_ndof;   ///< number of degrees of freedom
  ObjPropHandle<DOFMap> m_map; ///< mapping from indices to degrees of freedom

};

/** Handle for referring to a PhysObject. */
class PhysObjectHandle : public HandleBase
{
public:

  explicit PhysObjectHandle(int idx = -1) : HandleBase(idx) {}
};

#include "PhysObject.inl"

} // namespace BASim

#endif // PHYSOBJECT_HH

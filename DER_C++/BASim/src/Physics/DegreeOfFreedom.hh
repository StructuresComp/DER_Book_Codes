/**
 * \file DegreeOfFreedom.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/17/2009
 */

#ifndef DEGREEOFFREEDOM_HH
#define DEGREEOFFREEDOM_HH

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Core/Handle.hh"

namespace BASim {

/** Handle for referring to a degree of freedom.  */
class DofHandle : public HandleBase
{
public:

  /// Degrees of freedom can be associated to vertices, edges
  enum Type { VERTEX_DOF, EDGE_DOF };

  explicit DofHandle(int idx = -1)
    : HandleBase(idx)
    , m_dofNum(-1)
  {}

  DofHandle(const DofHandle& handle)
    : HandleBase(handle.m_idx)
    , m_type(handle.m_type)
    , m_handle(handle.m_handle)
    , m_dofNum(handle.m_dofNum)
  {}

  DofHandle& operator= (const DofHandle& handle)
  {
    m_idx = handle.m_idx;
    m_type = handle.m_type;
    m_handle = handle.m_handle;
    m_dofNum = handle.m_dofNum;
    return *this;
  }

  Type getType() const { return m_type; }
  void setType(const Type& type) { m_type = type; }

  int getNum() const { return m_dofNum; }
  void setNum(int num) { m_dofNum = num; }

  const HandleBase& getHandle() const { return m_handle; }
  void setHandle(const HandleBase& handle) { m_handle = handle; }

protected:

  Type m_type;
  HandleBase m_handle;
  int m_dofNum; ///< which degree of freedom
};

/** Mapping from degrees of freedom to indices. */
class DOFMap
{
public:

  DOFMap() {}
  ~DOFMap() {}

  void addMapping(const DofHandle& dof, int index)
  {
    m_dofToIndex.insert(std::make_pair(dof, index));
    m_indexToDof.insert(std::make_pair(index, dof));
  }

  const DofHandle& getDof(int index) const
  {
    std::map<int, DofHandle>::const_iterator it = m_indexToDof.find(index);
    assert(it != m_indexToDof.end());
    return it->second;
  }

  int getIndex(const DofHandle& handle) const
  {
    std::map<DofHandle, int>::const_iterator it = m_dofToIndex.find(handle);
    assert(it != m_dofToIndex.end());
    return it->second;
  }

protected:

  std::map<DofHandle, int> m_dofToIndex;
  std::map<int, DofHandle> m_indexToDof;
};

} // namespace BASim

#endif // DEGREEOFFREEDOM_HH

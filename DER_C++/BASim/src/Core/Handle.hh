/**
 * \file Handle.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef HANDLE_HH
#define HANDLE_HH

namespace BASim {

/** Base class for handles. */
class HandleBase
{
public:

  explicit HandleBase(int idx = -1) : m_idx(idx) {}

  virtual ~HandleBase() {}

  bool operator== (const HandleBase& rhs) const { return m_idx == rhs.m_idx; }
  bool operator!= (const HandleBase& rhs) const { return m_idx != rhs.m_idx; }
  bool operator< (const HandleBase& rhs) const { return m_idx < rhs.m_idx; }
  bool operator> (const HandleBase& rhs) const { return m_idx > rhs.m_idx; }

  int idx() const { return m_idx; }

  bool isValid() const { return m_idx != -1; }

  void __increment() { ++m_idx; }
  void __decrement() { --m_idx; }

protected:

  int m_idx;
};

} // namespace BASim

#endif // HANDLE_HH

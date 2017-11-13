#ifndef MATRIXBASE_HH
#define MATRIXBASE_HH

/**
 * \file MatrixBase.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/05/2009
 */

#include "BASim/src/Core/Definitions.hh"

namespace BASim {

/** Interface for matrices. */
class MatrixBase
{
public:

  enum Flags {
    SYMMETRIC         = 1,
    BANDED            = 1 << 1,
    TRIDIAGONAL       = 1 << 2,
    POSITIVE_DEFINITE = 1 << 3
  };

  virtual ~MatrixBase() {}

  /** Called before assigning matrix values. */
  virtual int begin() { return 0; }

  virtual Scalar operator() (int i, int j) const = 0;
  virtual int set(int i, int j, Scalar val) = 0;
  virtual int add(int i, int j, Scalar val) = 0;
  virtual int add(const IntArray& rowIdx, const IntArray& colIdx,
                  const MatXd& values) = 0;
  virtual int add(const IndexArray& rowIdx, const IndexArray& colIdx,
                  const MatXd& values) = 0;
  virtual int scale(Scalar val) = 0;
  virtual int setZero() = 0;
  virtual int zeroRows(const IntArray& idx, Scalar diag = 1.0) = 0;
  virtual int multiply(VecXd& y, Scalar s, const VecXd& x) const = 0;

  /**
   * Called after all matrix values have been set.
   */
  virtual int finalizeNonzeros() { return 0; }

  /**
   * Called if a new set of matrix values is going to be set.
   */
  virtual int resetNonzeros() { return 0; }
   
  /**
   * Some petsc stuff here
   */ 
  virtual int finalize() { return 0; }

  inline int rows() const { return m_rows; }
  inline int cols() const { return m_cols; }

  inline bool isFlagSet(const Flags& flag) const { return (m_flags & flag); }
  inline void setFlag(const Flags& flag) { m_flags = (Flags) (m_flags | flag); }
  inline void unsetFlag(const Flags& flag) { m_flags = (Flags) (m_flags & (~flag)); }

  virtual void print() const {}

protected:

  MatrixBase(int r, int c)
    : m_rows(r)
    , m_cols(c)
  {
    m_flags = (Flags) 0;
  }

  int m_rows;
  int m_cols;
  Flags m_flags;
};

} // namespace BASim

#endif // MATRIXBASE_HH

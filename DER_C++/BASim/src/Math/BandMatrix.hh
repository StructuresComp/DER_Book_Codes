#ifndef BANDMATRIX_HH
#define BANDMATRIX_HH

/**
 * \file BandMatrix.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 11/16/2009
 */

#include "BASim/src/Math/MatrixBase.hh"

namespace BASim {

/** Class for storing band matrices in a two-dimensional array. */
class BandMatrix : public MatrixBase
{
public:

  /** Creates an m-by-n band matrix with kl non-zero sub-diagonals and
      ku non-zero super-diagonals. */
  BandMatrix(int m, int n, int kl, int ku)
    : MatrixBase(m, n)
    , m_kl(kl)
    , m_ku(ku)
    , m_size((kl + ku + 1) * n)
  {
    m_data = new Scalar[m_size];
    setZero();

    m_lower = new int[m_rows];
    m_upper = new int[m_rows];
    for (int i = 0; i < m_rows; ++i) {
      m_lower[i] = std::max(i - m_kl, 0);
      m_upper[i] = std::min(i + m_ku + 1, m_cols);
    }
  }

  virtual ~BandMatrix()
  {
    delete [] m_data;
    delete [] m_lower;
    delete [] m_upper;
  }

  virtual Scalar operator() (int i, int j) const
  {
    if (!indicesValid(i, j)) return 0;
    return m_data[(m_ku + i - j) * cols() + j];
  }

  Scalar& operator() (int i, int j)
  {
    assert(indicesValid(i, j));
    return m_data[(m_ku + i - j) * cols() + j];
  }

  virtual int set(int i, int j, Scalar val)
  {
    (*this)(i, j) = val;
    return 0;
  }

  virtual int add(int i, int j, Scalar val)
  {
    (*this)(i, j) += val;
    return 0;
  }

  virtual int add(const IntArray& rowIdx, const IntArray& colIdx,
                  const MatXd& values)
  {
    size_t nr = rowIdx.size();
    size_t nc = colIdx.size();

    for (size_t i = 0; i < nr; ++i) {
      int r = rowIdx[i];
      assert(r >= 0);
      assert(r < rows());
      Scalar* val = &m_data[(m_ku + r) * m_cols];
      for (size_t j = 0; j < nc; ++j) {
        assert(indicesValid(r, colIdx[j]));
        *(val + colIdx[j] * (1 - m_cols)) += values(i, j);
      }
    }

    return 0;
  }

  virtual int add(const IndexArray& rowIdx, const IndexArray& colIdx,
                  const MatXd& values)
  {
    int nr = rowIdx.size();
    int nc = colIdx.size();

    for (int i = 0; i < nr; ++i) {
      int r = rowIdx[i];
      Scalar* val = &m_data[(m_ku + r) * m_cols];
      for (int j = 0; j < nc; ++j) {
        assert(indicesValid(r, colIdx[j]));
        *(val + colIdx[j] * (1 - m_cols)) += values(i, j);
      }
    }

    return 0;
  }

  virtual int scale(Scalar val)
  {
    for (int i = 0; i < m_size; ++i) m_data[i] *= val;
    return 0;
  }

  virtual int setZero()
  {
    for (int i = 0; i < m_size; ++i) m_data[i] = 0;
    return 0;
  }

  virtual int zeroRows(const IntArray& idx, Scalar diag = 1.0)
  {
    for (int i = 0; i < (int) idx.size(); ++i) {
      int r = idx[i];
      assert(r >=0 && r < m_rows);
      int lower = m_lower[r];
      int upper = m_upper[r];
      for (int j = lower; j < upper; ++j) {
        (*this)(r, j) = 0;
      }
      (*this)(r, r) = diag;
    }
    return 0;
  }

  virtual int multiply(VecXd& y, Scalar s, const VecXd& x) const
  {
    assert(y.size() == m_rows);
    assert(x.size() == m_cols);

    for (int i = 0; i < m_rows; ++i) {
      int lower = m_lower[i];
      int upper = m_upper[i];
      const Scalar* val = &m_data[(m_ku + i - lower) * m_cols + lower];
      Scalar sum = 0;
      for (int j = lower; j < upper; ++j, val += (1 - m_cols)) {
        sum += (*val) * x[j];
      }
      y[i] += s * sum;
    }
    return 0;
  }

  virtual void print() const
  {
    std::cout << "[";
    for (int i = 0; i < m_rows; ++i) {
      for (int j = 0; j < m_cols; ++j) {
        std::cout << (*this)(i, j);
        if (j < m_cols - 1) std::cout << ", ";
      }
      if (i < m_rows - 1) std::cout << "; ";
    }
    std::cout << "]" << std::endl;
  }

  /** Number of sub-diagonals */
  int kl() const { return m_kl; }

  /** Number of super-diagonals */
  int ku() const { return m_ku; }

  /** Returns underlying array of values */
  Scalar* data() { return m_data; }

protected:

  bool indicesValid(int r, int c) const
  {
    return ((r >= 0) && (r < rows()) && (c >= 0) && (c < cols())
            && (r - c <= m_kl) && (c - r <= m_ku));
  }

  int m_kl; ///< Number of sub-diagonals
  int m_ku; ///< Number of super-diagonals
  int m_size; ///< Size of the array holding the non-zero entries
  Scalar* m_data; ///< Array that stores the non-zero entries
  int* m_lower; ///< For each row, stores the smallest valid column index
  int* m_upper; ///< For each row, stores the largest valid column index
};

} // namespace BASim

#endif // BANDMATRIX_HH

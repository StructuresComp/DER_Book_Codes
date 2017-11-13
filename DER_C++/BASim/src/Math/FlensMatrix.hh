#ifndef FLENSMATRIX_HH
#define FLENSMATRIX_HH

/**
 * \file FlensMatrix.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/06/2009
 */

namespace BASim {

typedef flens::SparseGeMatrix< flens::CRS<Scalar> > SGeMatrix;

/** Wrapper around an FLENS matrix. */
template <class FM>
class FlensMatrix : public MatrixBase
{
public:

  FlensMatrix(int numRows, int numCols, int k=1)
    : m_matrix(numRows, numCols, k)
  {}

  Scalar& operator() (int i, int j)
  {
    return m_matrix(i + 1, j + 1);
  }

  const Scalar& operator() (int i, int j) const
  {
    return const_cast<FlensMatrix<FM>*>(this)->m_matrix(i + 1, j + 1);
  }

  void finalize()
  {
    m_matrix.finalize();
  }

  const FM& getFlensMatrix() const
  {
    return m_matrix;
  }

  FM& getFlensMatrix()
  {
    return m_matrix;
  }

protected:

  FM m_matrix;
};

} // namespace BASim

#endif // FLENSMATRIX_HH

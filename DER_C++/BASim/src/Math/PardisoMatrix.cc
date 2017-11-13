#include "PardisoMatrix.hh"

namespace BASim 
{

// NOTE ON IMPLEMENTATION
//
//  It is important to realize that with map, using the accessor like so
//
//  double x = matrix[i][j]
//
//  will automatically create a new entry in the map even if one does not exist yet.
//  This is true even if you are accessing (not writing) the map.
//
//  The correct way to iterate over elements in the map is to use iterators.
//  If you need to know whether there is a value mapped to a key, use the "find" function,
//  which will return the "end()" iterator if the value doesn't exist:
//
//  // row and cell are iterators
//  if ((row = matrix.find(i)) != matrix.end() {
//	  if ((cell = row.find(j)) != row.end()) {
//      ...
//    }
//  }

PardisoMatrix::PardisoMatrix( int rows, int cols )
: MatrixBase(rows,cols)
, m_finalized(false)
, m_nnz(-1)
, m_vals()
, m_colindices()
, m_rowstarts()
, m_matrix()
{
  assert( m_rows >= 0 );
  assert( m_cols >= 0 );
}

PardisoMatrix::~PardisoMatrix()
{
}

bool PardisoMatrix::isFinalized() const 
{ 
  return m_finalized; 
}

int PardisoMatrix::nonZeros() const 
{
  assert( m_finalized );
  
  return m_nnz; 
}

VecXd& PardisoMatrix::getVals()
{
  assert( m_finalized );
  assert( m_vals.size() == m_nnz );
  return m_vals;
}

Eigen::VectorXi& PardisoMatrix::getRowstarts() 
{
  assert( m_finalized );
  assert( m_rowstarts.size() == m_rows+1 );
  return m_rowstarts;
}

Eigen::VectorXi& PardisoMatrix::getColindices() 
{
  assert( m_finalized );
  assert( m_colindices.size() == m_vals.size() );
  assert( m_colindices.size() == m_nnz );
  return m_colindices;
}

//PardisoMatrix* PardisoMatrix::getTranspose() const
//{
//  assert( !m_finalized );
//
//  // Create a new Pardiso Matrix with the same dimensions as this matrix.
//  PardisoMatrix* pm = new PardisoMatrix(m_cols,m_rows);
//  //pm->m_nnz = m_nnz;
//  
//  // For each row of this matrix
//  for( PMatRows::const_iterator i = m_matrix.begin(); i != m_matrix.end(); ++i )
//  {
//    PMatCols cols = (*i).second;
//    // For each column of this matrix
//    for( PMatCols::const_iterator j = cols.begin(); j != cols.end(); ++j ) 
//    {
//      int rownum = (*i).first;
//      assert( rownum >= 0 ); assert( rownum < m_rows );
//      int colnum = (*j).first;
//      assert( colnum >= 0 ); assert( colnum < m_cols );
//      double val = (*j).second;
//      pm->m_matrix[colnum][rownum] = val;
//    }
//  }
//
//  return pm;
//}

// TODO: Rearrange the computation here, avoid allocating a vector and then copying to an array
int PardisoMatrix::finalizeNonzeros()
{
  if( !m_finalized )
  {
    assert( m_vals.size() == 0 );
    assert( m_colindices.size() == 0 );
    assert( m_rowstarts.size() == 0 );
    
    std::vector<double> v_vals;
    std::vector<int> v_colindices;
    std::vector<int> v_rowstarts;
    
    //int valindex = 0;
    //rowvectors = new vector<pair<int,int> >[m_rows];
    
    // For each row
    for( std::map<int, std::map<int, double> >::iterator i = m_matrix.begin(); i != m_matrix.end(); ++i )
    {
      std::map<int, double> cols = (*i).second;
      bool firsttime = true;
      // For each column
      for( std::map<int, double>::iterator j = cols.begin(); j != cols.end(); ++j )
      {
        int rownum = (*i).first;
        assert( rownum >= 0 ); assert( rownum < m_rows );
        int colnum = (*j).first;
        assert( colnum >= 0 ); assert( colnum < m_cols );
        double val = (*j).second;
        if( firsttime )
        {
          firsttime = false;
          v_rowstarts.push_back(v_colindices.size());
        }
        v_vals.push_back(val);
        v_colindices.push_back(colnum);
        
        // prepare hashmap for addScaled
        // this is a macro!
        //all_index_hashes.resize(_rows);
        //all_index_hashes[rownum][colnum] = valindex;
        
        // prepare array of row vectors for multVector
        //pair<int, int> rowvectval(colnum, valindex);
        //rowvectors[rownum].push_back(rowvectval);
        
        // increment valindex
        //++valindex;
      }
    }
    v_rowstarts.push_back(v_colindices.size());
    
    m_nnz = v_vals.size();
    
    assert( v_vals.size() == v_colindices.size() );
    assert( (int) v_rowstarts.size() == m_rows+1 );
    assert( v_rowstarts[m_rows] == m_nnz );
    
    /*
     // should I resize?
     for (int i = 0 ; i < _rows ; i++) {
     all_index_hashes[i].resize(40);
     }
     */
    
    // If m_nnz == 0, the matrix is clear -- just don't allocate anything.
    if( m_nnz > 0 ) 
    {
      m_vals.resize(v_vals.size());
      m_colindices.resize(v_colindices.size());
      m_rowstarts.resize(v_rowstarts.size());
      
      for( unsigned int i = 0; i < v_vals.size(); ++i ) 
      {
        m_vals(i) = v_vals[i];
        m_colindices(i) = v_colindices[i];
      }
      for( unsigned int i = 0; i < v_rowstarts.size(); ++i ) 
      {
        m_rowstarts(i) = v_rowstarts[i];
      }
    }
    m_matrix.clear();
    m_finalized = true;
  }
  
  return 0;
}

int PardisoMatrix::resetNonzeros()
{
  m_finalized = false;
  m_nnz = -1;
  m_vals.resize(0);
  m_colindices.resize(0);
  m_rowstarts.resize(0);
  m_matrix.clear();
  
  return 0;
}
  
double PardisoMatrix::operator()( int i, int j ) const
{
  assert( i >= 0 ); assert( i < m_rows );
  assert( j >= 0 ); assert( j < m_cols );
  
  if( !m_finalized )
  {
    // If the requested row is empty, value hasn't been set
    std::map<int, std::map<int, double> >::const_iterator row = m_matrix.find(i);
    if( row == m_matrix.end() ) return 0.0;
    
    // If request col does not exit, value hasn't been set
    std::map<int, double>::const_iterator col = (*row).second.find(j);
    if( col == (*row).second.end() ) return 0.0;
    
    // Otherwise we have it!
    return (*col).second;
  }
  
  int valindex = findValueIndex(i,j);
  assert( valindex >= -1 ); assert( valindex < m_nnz );
  if( valindex < 0 ) return 0.0;
  return m_vals(valindex);
}

int PardisoMatrix::set( int i, int j, double val )
{
  assert( !m_finalized );
  assert( i >= 0 ); assert( i < m_rows );
  assert( j >= 0 ); assert( j < m_cols );
  
  m_matrix[i][j] = val;
  
  return 0;
}

int PardisoMatrix::add( int i, int j, double val )
{
  assert( i >= 0 ); assert( i < m_rows );
  assert( j >= 0 ); assert( j < m_cols );
  
  if( !m_finalized )
  {
    m_matrix[i][j] += val;
  } 
  else
  {
    int valindex = findValueIndex(i,j);
    assert( valindex >= -1 ); assert( valindex < m_nnz );
    
    // Requested element exists in finalized non-zero structure
    if( valindex != -1 )
    {
      m_vals(valindex) += val;
    }
    // Requested element does not exist in finalized non-zero structure
    else
    {
      std::cerr << "\033[31;1mERROR IN PARDISOMATRIX:\033[m attempted to access non-existant element in finalized matrix using add( int i, int j, double val ). Exiting." << std::endl;
      exit(0);
    }
  }
  
  return 0;
}

int PardisoMatrix::add( const IntArray& rowIdx, const IntArray& colIdx, const MatXd& values )
{
  int nr = rowIdx.size();
  int nc = colIdx.size();
  
  assert( values.rows() == nr );
  assert( values.cols() == nc );
  
  if( !m_finalized )
  {  
    for( int i = 0; i < nr; ++i ) for( int j = 0; j < nc; ++j ) 
    {
      assert( rowIdx[i] >= 0 ); assert( rowIdx[i] < rows() );
      assert( colIdx[j] >= 0 ); assert( colIdx[j] < cols() );
      m_matrix[rowIdx[i]][colIdx[j]] += values(i,j);
    }
  }
  else
  {    
    for( int i = 0; i < nr; ++i ) for( int j = 0; j < nc; ++j ) 
    {
      assert( rowIdx[i] >= 0 ); assert( rowIdx[i] < rows() );
      assert( colIdx[j] >= 0 ); assert( colIdx[j] < cols() );
      int valindex = findValueIndex(rowIdx[i],colIdx[j]);
      assert( valindex >= -1 ); assert( valindex < m_nnz );
      if( valindex < 0 )
      {
        std::cerr << "\033[31;1mERROR IN PARDISOMATRIX:\033[m attempted to access non-existant element in finalized matrix using add( const IntArray& rowIdx, const IntArray& colIdx, const MatXd& values ). Exiting." << std::endl;
        exit(0);
      }
      m_vals(valindex) += values(i,j);
    }
  }
  
  return 0;
}

int PardisoMatrix::add( const IndexArray& rowIdx, const IndexArray& colIdx, const MatXd& values )
{
  int nr = rowIdx.size();
  int nc = colIdx.size();
  
  assert( values.rows() == nr );
  assert( values.cols() == nc );
  
  if( !m_finalized )
  {  
    for( int i = 0; i < nr; ++i ) for( int j = 0; j < nc; ++j ) 
    {
      assert( rowIdx[i] >= 0 ); assert( rowIdx[i] < rows() );
      assert( colIdx[j] >= 0 ); assert( colIdx[j] < cols() );
      m_matrix[rowIdx[i]][colIdx[j]] += values(i,j);
    }
  }
  else
  {    
    for( int i = 0; i < nr; ++i ) for( int j = 0; j < nc; ++j ) 
    {
      assert( rowIdx[i] >= 0 ); assert( rowIdx[i] < rows() );
      assert( colIdx[j] >= 0 ); assert( colIdx[j] < cols() );
      int valindex = findValueIndex(rowIdx[i],colIdx[j]);
      assert( valindex >= -1 ); assert( valindex < m_nnz );
      if( valindex < 0 )
      {
        std::cerr << "\033[31;1mERROR IN PARDISOMATRIX:\033[m attempted to access non-existant element in finalized matrix using add( const IndexArray& rowIdx, const IndexArray& colIdx, const MatXd& values ). Exiting." << std::endl;
        exit(0);
      }
      m_vals(valindex) += values(i,j);
    }
  }
  
  return 0;  
}

int PardisoMatrix::scale( double val )
{
  if( !m_finalized )
  {
    std::map<int, std::map<int, double> >::iterator rowitr;
    for( rowitr = m_matrix.begin(); rowitr != m_matrix.end(); ++rowitr )
    {
      std::map<int, double>::iterator colitr;
      for( colitr = (*rowitr).second.begin(); colitr != (*rowitr).second.end(); ++colitr )
      {
        (*colitr).second *= val;
      }
    }
  }
  else 
  {
    m_vals *= val;
  }
  
  return 0;  
}

int PardisoMatrix::setZero()
{
  return scale(0.0);
}

int PardisoMatrix::zeroRows( const IntArray& idx, double diag )
{
#ifdef DEBUG
  // Ensure the user provided valid rows
  for( int i = 0; i < (int) idx.size(); ++i ) 
  {
    assert( idx[i] >= 0 );
    assert( idx[i] < rows() );
  }
#endif
  
  if( !m_finalized ) 
  {
    for( int i = 0; i < (int) idx.size(); ++i )
    {
      // If the requested row has any elements
      std::map<int, std::map<int, double> >::iterator rowitr;
      rowitr = m_matrix.find( idx[i] );
      if( rowitr == m_matrix.end() ) continue;
      
      // Zero all but the diagonal of that row
      std::map<int, double>::iterator colitr;
      for( colitr = (*rowitr).second.begin(); colitr != (*rowitr).second.end(); ++colitr )
      {
        // If off diagonal
        if( (*colitr).first != idx[i] ) 
        {
          (*colitr).second = 0.0;
        }
        // Else, if on the diagonal
        else
        {
          (*colitr).second = diag;
        }
      }
    }
  } 
  else
  {
    assert( m_rowstarts.size() == rows()+1 );
    
    // For each row specified in idx
    for( int i = 0; i < (int) idx.size(); ++i ) 
    {
      int rowstart = m_rowstarts[idx[i]];
      int onepastrowend = m_rowstarts[idx[i]+1];
      assert( rowstart <= onepastrowend );
      
      for( int j = rowstart; j < onepastrowend; ++j )
      {
        if( m_colindices[j] != idx[i] )
        {
          m_vals(j) = 0.0;
        }
        else 
        {
          m_vals(j) = diag;
        }
      }
    }
  }
  
  return 0;
}

int PardisoMatrix::multiply( VecXd& y, double s, const VecXd& x ) const
{
  assert( y.size() == x.size() );
  assert( y.size() == rows() );
  
  if( m_finalized ) 
  {
    // For each row of this matrix
    for( int row = 0; row < rows(); ++row )
    {
      // Compute dot product of the current row with the input vector
      double sum = 0.0;
      // For each column in the current row
      for( int j = m_rowstarts[row]; j < m_rowstarts[row+1]; ++j )
      {
        int col = m_colindices[j];
        assert( col >= 0 ); assert( col < cols() );
        sum += m_vals[j]*x[col];
      }
      y[row] += s*sum;
    }
  }
  else 
  {
    // For each row
    for( std::map<int, std::map<int, double> >::const_iterator i = m_matrix.begin(); i != m_matrix.end(); ++i )
    {
      std::map<int, double> rowcols = (*i).second;
      int rownum = (*i).first;
      assert( rownum >= 0 ); assert( rownum < rows() );
      double total = 0.0;
      // For each column
      for( std::map<int, double>::const_iterator j = rowcols.begin(); j != rowcols.end(); ++j )
      {
        int colnum = (*j).first;
        assert( colnum >= 0 ); assert( colnum < cols() );
        double val = (*j).second;
        total += val * x[colnum];
      }
      y(rownum) += s*total;
    }
  }
  
  return 0;  
}

void PardisoMatrix::print()
{
  if( !m_finalized )
  {
    // For each row
    for( int i = 0; i < rows(); ++i )
    {        
      std::map<int, std::map<int, double> >::const_iterator row = m_matrix.find(i);
      // If the row is empty
      if( row == m_matrix.end() )
      {
        for( int j = 0; j < cols(); ++j ) std::cout << "X   ";
      }
      // Else if the row is not empty
      else
      {
        for( int j = 0; j < cols(); ++j )
        {
          std::map<int, double>::const_iterator col = (*row).second.find(j);
          if( col == (*row).second.end() ) std::cout << "X   ";
          else std::cout << (*col).second << "   ";        
        }
      }
      std::cout << std::endl;
    }
  }
  else
  {
    // For each row
    for( int i = 0; i < rows(); ++i )
    {
      int rowstart = m_rowstarts[i];
      int onepastrowend = m_rowstarts[i+1];
      for( int j = 0; j < cols(); ++j )
      {
        if( (rowstart<onepastrowend)&&(m_colindices[rowstart]==j) )
        {
          std::cout << m_vals[rowstart] << "   ";
          ++rowstart;
        }
        else
        {
          std::cout << "X   ";
        }
      }
      std::cout << std::endl;
    }
  }
}
  
void PardisoMatrix::runSanityChecks()
{
  if( m_finalized )
  {
    assert( m_nnz == m_vals.size() );
    assert( m_nnz == m_colindices.size() );
    assert( m_rows == m_rowstarts.size()-1 );
    for( int i = 0; i < m_colindices.size(); ++i ) 
    {
      assert( m_colindices[i] >= 0 );
      assert( m_colindices[i] < m_cols );
    }
    for( int i = 0; i < m_rowstarts.size()-1; ++i )
    {
      assert( m_rowstarts[i] <= m_rowstarts[i+1] );
    }
    assert( m_rowstarts[m_rowstarts.size()-1] == m_nnz );
  }
  else
  {
  }
}  

int PardisoMatrix::findValueIndex( const int& row, const int& col ) const
{
  assert( row >= 0 ); assert( row < m_rows );
  assert( col >= 0 ); assert( col < m_cols );
  
  assert( m_finalized );
  
  int start_of_row = m_rowstarts(row);
  assert( start_of_row >= 0 );
  //assert( start_of_row < m_nnz );
  
  int start_of_next_row = m_rowstarts(row+1);
  assert( start_of_next_row >= 0 );
  //assert( start_of_next_row < m_nnz );
  
  for( int i = start_of_row; i < start_of_next_row; ++i )
  {
    assert( m_colindices(i) >= 0 );
    assert( m_colindices(i) < m_cols );
    if( m_colindices(i) == col ) 
    {
      // Found the element
      return i;
    }
  }
  
  // Failed to find the element
  return -1;
}

}

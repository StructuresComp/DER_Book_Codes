#ifndef PARDISOMATRIX_HH
#define PARDISOMATRIX_HH

/**
 * \file PardisoMatrix.hh
 *
 * I think this was originally authored by either: R. Goldenthal, B. Audoly, or J. Sherrer :). 
 * Hacked up, cleaned up, and plugged into BASim by:
 * \author smith@cs.columbia.edu
 * \date 05/28/2010
 */

// TODO: Although this was written for use with Pardiso, its use is not restricted to pardiso alone.
//       Make a little more general, clean up a bit. 

#include "BASim/src/Math/MatrixBase.hh"

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cassert>
#include <cmath>

namespace BASim 
{
  class PardisoMatrix : public MatrixBase
  {
    
  public:

    /** 
     * Constructs a matrix of dimension rows x cols.
     *
     * @param rows Number of rows in the matrix
     * @param cols Number of cols in the matrix
     */
    PardisoMatrix( int rows, int cols );

    /** 
     * Destructor.
     */
    virtual ~PardisoMatrix();

    /**
     * True if the non-zero structure of the matrix has been finalized,
     * false otherwise.
     */
    bool isFinalized() const;
    
    /**
     * Returns the number of nonzeros in the finalized matrix.
     * Do not call before finalizing the structure. 
     */
    int nonZeros() const;
    
    /**
     * Returns a pointer to an array of the matrix's contents suitable for use
     * with Pardiso. 
     */
    VecXd& getVals();
    
    /**
     * Returns a pointer to an array containing the indices of the first element of each row.
     * That is, getRowStarts()[i] points to the first element of row i in getVals().
     */
    Eigen::VectorXi& getRowstarts();
    
    /**
     * Returns a pointer to an array containing the column index of each element.
     * That is, getColindices()[i] is the column of the ith element of getVals(). 
     */
    Eigen::VectorXi& getColindices();
    
    /**
     * Finalizes the non-zero structure of the matrix, prepare for use by Pardiso. 
     */
    int finalizeNonzeros();
    
    /**
     * Returns to the matrix to a completely blank slate as if constructor had been called. 
     */
    int resetNonzeros();
    
    /**
     * Returns a value as if this were a dense matrix; if A(i,j) has been set, returns
     * the value of the element. If A(i,j) hasn't been set, returns 0. One CAN NOT set 
     * with this method.
     *
     * @param i ith row of the matrix. 
     * @param j jth column of the matrix.
     */
    double operator()( int i, int j ) const;
    
    /**
     * Inserts scalar value a into the specfied row and column.
     * Do not call after calling finalize.
     *
     * @param i ith row of the matrix. 
     * @param j jth column of the matrix.
     * @param val Value to insert at position i,j. 
     */
    int set( int i, int j, double val );
    
    /**
     * If value i,j was already set, adds val to that element. 
     * If value i,j was not set, sets i,j to val. 
     *
     * @param i ith row of the matrix. 
     * @param j jth column of the matrix.
     * @param val Value to add to position i,j. 
     */
    int add( int i, int j, double val );
    
    /**
     * Adds values in matrix values to this matrix. values(i,j) is added to
     * position rowIdx[i] and colIdx[i] in this matrix.
     *
     * @param rowIdx ith row of values is rowIdx(i)th row of this matrix.
     * @param colIdx jth column of values is colIdx(i)th row of this matrix.
     * @param values values to add to this matrix. 
     */
    int add( const IntArray& rowIdx, const IntArray& colIdx, const MatXd& values );
    
    /**
     * Adds values in matrix values to this matrix. values(i,j) is added to
     * position rowIdx[i] and colIdx[i] in this matrix.
     *
     * @param rowIdx ith row of values is rowIdx(i)th row of this matrix.
     * @param colIdx jth column of values is colIdx(i)th row of this matrix.
     * @param values values to add to this matrix. 
     */
    int add( const IndexArray& rowIdx, const IndexArray& colIdx, const MatXd& values );
    
    /**
     * Scales every value of this matrix by val.
     *
     * @param val value to scale each element of this matrix by.
     */  
    int scale( double val );
    
    /**
     * Sets all elements of the matrix to 0. Does not change the sparsity structure. 
     */  
    int setZero();
    
    /**
     * Zeros the off-diagonal elements of rows sepcified in idx, sets the diagonal
     * elements to diag. Does not modify the sparsity structure.
     *
     * @param idx array of row indices.
     * @param diag value to set diagonal elements to
     */    
    int zeroRows( const IntArray& idx, double diag = 1.0 );
    
    /**
     * Computes the matrix-vector multiplication y = s*this*x
     *
     * @param y vector to store s*this*x in
     * @param s value to scale the matrix-vector multiplication by
     * @param x vector to store result in
     */  
    int multiply( VecXd& y, double s, const VecXd& x ) const;
    
    void print();
    
    /**
     * Asserts some invariants.
     */
    void runSanityChecks();
    
  private:
    
    /**
     * Given a row and a column, returns the index into the vals array that refers
     * to the item, or -1 if it is not present.
     *
     *  @param row A row of the matrix
     *  @param col A column of the matrix
     */
    int findValueIndex( const int& row, const int& col ) const;
    
    // True if the matrix's non-zero structure has been finalized
    bool m_finalized;
    
    // Number of non-zero elements in the matrix
    int m_nnz;
    
    // Values of nonzero elements in row-major order
    VecXd m_vals;
    // Column indices of nonzer elements in row-major order
    Eigen::VectorXi m_colindices;
    // rowstarts[i] is the start of row i in colindices and vals
    Eigen::VectorXi m_rowstarts;
    
    std::map<int, std::map<int, double> > m_matrix;
  };  
}

#endif

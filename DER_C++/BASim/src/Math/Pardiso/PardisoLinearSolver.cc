/**
 * \file PardisoLinearSolver.cc
 *
 * \author smith@cs.columbia.edu
 * \date 05/29/2010
 */

#include "PardisoLinearSolver.hh"

using namespace std;

namespace BASim {

PardisoLinearSolver::PardisoLinearSolver( PardisoMatrix& A )
: LinearSolverBase(A)
, prdsomat(A)
{
  assert(A.rows() == A.cols());
}

PardisoLinearSolver::~PardisoLinearSolver()
{
}

void PardisoLinearSolver::parsePardisoError( int error ) const 
{
  switch( error )
  {
    case 0:
    {
      //std::cout << "No error." << std::endl;
      break;
    }
    case -1:
    {
      std::cerr << "Input inconsistent. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -2:
    {
      std::cerr << "Not enough memory. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -3:
    {
      std::cerr << "Reordering problem. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -4:
    {
      std::cerr << "Zero pivot, numerical fact. or iterative refinement problem. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -5:
    {
      std::cerr << "Unclassified (internal) error. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -6:
    {
      std::cerr << "Preordering failed (matrix types 11, 13 only). Exiting." << std::endl;
      exit(0);
      break;
    }
    case -7:
    {
      std::cerr << "Diagonal matrix problem. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -8:
    {
      std::cerr << "32-bit integer overflow problem. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -10:
    {
      std::cerr << "No license file pardiso.lic found. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -11:
    {
      std::cerr << "License is expired. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -12:
    {
      std::cerr << "Wrong username or hostname. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -100:
    {
      std::cerr << "Reached maximum number of Krylov-subspace iteration in iterative solver. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -101:
    {
      std::cerr << "No sufficient convergence in Krylov-subspace iteration within 25 iterations. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -102:
    {
      std::cerr << "Error in Krylov-subspace iteration. Exiting." << std::endl;
      exit(0);
      break;
    }
    case -103:
    {
      std::cerr << "Break-Down in Krylov-subspace iteration. Exiting." << std::endl;
      exit(0);
      break;
    }
    default:
    {
      std::cerr << "INVALID ERROR CODE ENCOUNTERED. THIS IS A BUG. EXITING." << std::endl;
      exit(0);
      break;
    }
  }
}


int PardisoLinearSolver::solve( VecXd& x, const VecXd& b )
{
  assert( prdsomat.isFinalized() );
  #ifdef DEBUG
    prdsomat.runSanityChecks();
  #endif
  
  // 1: Phase 1: Fill-reduction analysis and symbolic factorization
  // 2: Phase 2: Numerical factorization
  // 3: Phase 3: Forward and Backward solve including iterative refinement
  // 4: Termination and Memory Release Phase (Phase <= 0)
  
  assert( prdsomat.rows() == prdsomat.cols() );

  // Extract necessary data from the sparse matrix
  //int nnz = sparseMatrix.nonZeros();
  Eigen::VectorXi& rowstarts = prdsomat.getRowstarts();
  Eigen::VectorXi& colindices = prdsomat.getColindices();
  VecXd& vals = prdsomat.getVals();

  int n = prdsomat.rows();

  // Real and nonsymmetric
  int mtype = 11;

  // Number of right hand sides
  int nrhs = 1;

  // Internal solver memory pointer
  void* pt[64];

  // Pardiso control parameters
  int iparm[64];
  double dparm[64];

  // Double dummy
  double ddum;
  // Integer dummy
  int idum;

  // If 0, use sparse direct solver. If 1, use multi-recursive iterative solver. 
  int solver = 0;

  /////////////////////////////////////////////////////////////////////////////
  // ENSURE THE LICENSE FILE IS PRESENT

  // Track errors
  int error = 0;
  F77_FUNC(pardisoinit)( pt,  &mtype, &solver, iparm, dparm, &error );
  parsePardisoError( error );


  /////////////////////////////////////////////////////////////////////////////
  // REORDERING AND SYMBOLIC FACTORIZATION

  // Determine the desired number of threads to use.
  int num_threads = 1; // omp_get_max_threads(); // fromString<int>(getenv("OMP_NUM_THREADS"));
  //std::cout << "NUM_THREADS: " << num_threads << std::endl;
  iparm[2] = num_threads;

  // Number of factors with identical nonzero sparsity structure that the user would like to keep in memory.
  int maxfct = 1;
  //Actual matrix for the solution phase (matrix the user would like to factorize). 1 <= mnum <= maxfct
  int mnum = 1;
  // If 0, pardiso generates no output. If 1, pardiso prints statistical info.
  int msglvl = 0;

  // Convert matrix from 0-based C-notation to Fortran 1-based notation.
  rowstarts.cwise() += 1;
  colindices.cwise() += 1;

  // Reordering and Symbolic Factorization.  This step also allocates all memory that is necessary for the factorization.
  int phase = 11;
  error = 0;
  F77_FUNC(pardiso)( pt, &maxfct, &mnum, &mtype, &phase, &n, vals.data(), rowstarts.data(), colindices.data(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm );
  parsePardisoError( error );

  //std::cout << "  Reordering completed." << std::endl;
  //std::cout << "    Number of nonzeros in factors = " << iparm[17] << std::endl;
  //std::cout << "    Number of factorization MFLOPS = " << iparm[18] << std::endl;

  // Numerical factorization
  phase = 22;
  F77_FUNC(pardiso)( pt, &maxfct, &mnum, &mtype, &phase, &n, vals.data(), rowstarts.data(), colindices.data(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm );
  parsePardisoError( error );

  //std::cout << "  Factorization completed." << std::endl;


  /////////////////////////////////////////////////////////////////////////////
  // BACK SUBSTITUTION AND ITERATIVE REFINEMENT

  // Back substitution and iterative refinement
  phase = 33;
  // Max numbers of iterative refinement steps
  iparm[7] = 1;
  VecXd nonconstrhs = const_cast<VecXd&>(b);
  F77_FUNC(pardiso)( pt, &maxfct, &mnum, &mtype, &phase, &n, vals.data(), rowstarts.data(), colindices.data(), &idum, &nrhs, iparm, &msglvl, nonconstrhs.data(), x.data(), &error, dparm );
  parsePardisoError( error );

  // Convert matrix back to 0-based C-notation.
  rowstarts.cwise() -= 1;
  colindices.cwise() -= 1;

  // Release internal memory
  phase = -1;

  F77_FUNC(pardiso)( pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, rowstarts.data(), colindices.data(), &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm );
  parsePardisoError( error );

  return 0;
}

} // namespace BASim

/**
 * \file PardisoLinearSolver.hh
 *
 * \author smith@cs.columbia.edu
 * \date 05/29/2010
 */
 
// TODO: Add extra checks from Pardiso's return status (check the residual, etc)

#ifndef PARDISOLINEARSOLVER_HH
#define PARDISOLINEARSOLVER_HH

#include "BASim/src/Math/LinearSolverBase.hh"
#include "BASim/src/Math/PardisoMatrix.hh"

#define F77_FUNC(func)  func ## _

// PARDISO prototypes
extern "C" int F77_FUNC(pardisoinit)
(void *, int *, int *, int *, double *, int *);

extern "C" int F77_FUNC(pardiso)
(void *, int *, int *, int *, int *, int *, 
 double *, int *, int *, int *, int *, int *, 
 int *, double *, double *, int *, double *);

namespace BASim
{

/** 
 * Solves a sparse linear system using Pardiso. 
 */
class PardisoLinearSolver : public LinearSolverBase
{
public:

  PardisoLinearSolver( PardisoMatrix& A );
  virtual ~PardisoLinearSolver();

  void parsePardisoError( int error ) const;

  virtual int solve( VecXd& x, const VecXd& b );

private:
  PardisoMatrix& prdsomat;
};

} // namespace BASim

#endif // PARDISOLINEARSOLVER_HH

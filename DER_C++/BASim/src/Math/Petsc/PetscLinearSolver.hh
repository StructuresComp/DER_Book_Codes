/**
 * \file PetscLinearSolver.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/08/2009
 */

#ifndef PETSCLINEARSOLVER_HH
#define PETSCLINEARSOLVER_HH

#include <petscksp.h>
#include "BASim/src/Math/LinearSolverBase.hh"

namespace BASim {

/** Linear solver that uses PETSc. */
class PetscLinearSolver : public LinearSolverBase
{
public:

  PetscLinearSolver(MatrixBase& A)
    : LinearSolverBase(A)
    , m_x(NULL)
    , m_b(NULL)
  {
    assert(m_A.rows() == m_A.cols());

    // create solver context and associate system matrix with it
    KSPCreate(PETSC_COMM_SELF, &m_kspSolver);
    Mat& pA = smart_cast<PetscMatrix&>(m_A).getPetscMatrix();
    KSPSetOperators(m_kspSolver, pA, pA, SAME_NONZERO_PATTERN);
    KSPSetFromOptions(m_kspSolver);
    /*
    const KSPType kspType;
    KSPGetType(m_kspSolver, &kspType);
    std::cout << "KSP type is " << kspType << std::endl;

    PC pc;
    KSPGetPC(m_kspSolver, &pc);
    const PCType pcType;
    PCGetType(pc, &pcType);
    std::cout << "PC type is " << pcType << std::endl;
    */
    // create PETSc solution and right-hand-side vectors
    VecCreateSeq(PETSC_COMM_SELF, m_A.rows(), &m_x);
    VecCreateSeq(PETSC_COMM_SELF, m_A.rows(), &m_b);
  }

  ~PetscLinearSolver()
  {
    KSPDestroy(m_kspSolver);
    VecDestroy(m_x);
    VecDestroy(m_b);
  }

  int solve(VecXd& x, const VecXd& b)
  {
    PetscUtils::copyToPetscVector(m_b, b);
    PetscUtils::copyToPetscVector(m_x, x);

    KSPSolve(m_kspSolver, m_b, m_x);
    /*PetscInt its;
    KSPGetIterationNumber(m_kspSolver, &its);
    std::cout << "converged in " << its << " iterations" << std::endl;
    */

    PetscUtils::copyFromPetscVector(x, m_x);

    if (checkConvergedReason() == -1) return -1;
    return 0;
  }

  int checkConvergedReason()
  {
    KSPConvergedReason reason;
    KSPGetConvergedReason(m_kspSolver, &reason);

    if (reason == KSP_DIVERGED_NULL)
      std::cout << "Diverged null" << std::endl;
    else if (reason == KSP_DIVERGED_ITS)
      std::cout << "Diverged its" << std::endl;
    else if (reason == KSP_DIVERGED_NAN)
      std::cout << "Diverged nan" << std::endl;
    else if (reason == KSP_DIVERGED_BREAKDOWN_BICG)
      std::cout << "Diverged breakdown bicg" << std::endl;
    else if (reason == KSP_DIVERGED_DTOL)
      std::cout << "Diverged dtol" << std::endl;
    else if (reason == KSP_DIVERGED_INDEFINITE_PC)
      std::cout << "Diverged indefinite pc" << std::endl;
    else if (reason == KSP_DIVERGED_INDEFINITE_MAT)
      std::cout << "Diverged indefinite mat" << std::endl;
    else if (reason == KSP_DIVERGED_NONSYMMETRIC)
      std::cout << "Diverged nonsymmetric" << std::endl;
    else if (reason == KSP_DIVERGED_BREAKDOWN)
      std::cout << "Diverged breakdown" << std::endl;

    if (reason < 0) return -1;
    return 0;
  }

protected:

  KSP m_kspSolver;
  Vec m_x;
  Vec m_b;
};

} // namespace BASim

#endif // PETSCLINEARSOLVER_HH

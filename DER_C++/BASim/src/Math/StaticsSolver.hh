/**
 * \file StaticsSolver.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 05/14/2010
 */

#ifndef STATICSSOLVER_HH
#define STATICSSOLVER_HH

#include "BASim/src/Math/TimeSteppingBase.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/LinearSolverBase.hh"
#include "BASim/src/Math/SolverUtils.hh"
#include "BASim/src/Core/Timer.hh"

namespace BASim {

/** This class implements a statics solver. */
template <class ODE>
class StaticsSolver : public DiffEqSolver
{
public:

  StaticsSolver(ODE& ode)
    : m_diffEq(ode)
    , m_initial_residual(0)
    , m_residual(0)
    , m_curit(0)
  {
    m_A = m_diffEq.createMatrix();
    m_solver = SolverUtils::instance()->createLinearSolver(m_A);
  }

  ~StaticsSolver()
  {
    delete m_A;
    delete m_solver;
  }

  void execute()
  {
    m_diffEq.startStep();

    resize();
    setZero();
    m_diffEq.getScriptedDofs(m_fixed, m_desired);

    // copy initial guess for positions
    for (int i = 0; i < x0.size(); ++i) {
      x0(i) = m_diffEq.getX(i);
    }

    // computeResidual also sets m_rhs = F
    m_initial_residual = computeResidual();
    //std::cout << "initial residual: " << m_initial_residual << std::endl;
    m_increment[0] = 1;
    if ( isConverged() ) {
      m_diffEq.endStep();
      return;
    }
    m_increment.setZero();
    int m_curit = 0;

    for (; m_curit < m_maxit; ++m_curit) {
      START_TIMER("setup");
      // m_rhs = -F
      m_rhs *= -1.0;
      m_diffEq.startIteration();

      // m_A = dF/dx
      m_diffEq.evaluatePDotDX(1.0, *m_A);
      m_A->finalize();

      for (size_t i = 0; i < m_fixed.size(); ++i) {
        int idx = m_fixed[i];
        m_rhs[idx] = m_desired[i] - m_diffEq.getX(idx);
      }
      m_A->zeroRows(m_fixed);
      STOP_TIMER("setup");

      START_TIMER("solver");
      int status = m_solver->solve(m_increment, m_rhs);
      if( status < 0 ) std::cout << "\033[31;1mWARNING IN STATICSSOLVER:\033[m Problem during linear solve detected. " << std::endl;
      STOP_TIMER("solver");

      START_TIMER("setup");
      m_diffEq.increment_q(m_increment);
      x0 += m_increment;

      m_diffEq.endIteration();

      // check for convergence
      if ( isConverged() ) break;
      if (m_curit == m_maxit - 1) break;

      m_increment.setZero();
      m_A->setZero();
      STOP_TIMER("setup");
    }

    //std::cout << m_curit << std::endl;
    if (m_curit == m_maxit - 1)
    {
      std::cout << "\033[31;1mWARNING IN STATICSSOLVER:\033[m Newton solver failed to converge in max iterations: " << m_maxit << std::endl;
    }

    m_diffEq.endStep();
  }

  std::string getName() const
  {
    return "Statics solver";
  }

  void resize()
  {
    x0.resize(m_diffEq.ndof());
    m_rhs.resize(m_diffEq.ndof());
    m_lhs.resize(m_diffEq.ndof());
    m_increment.resize(m_diffEq.ndof());
    if (m_A->rows() != m_diffEq.ndof()) {
      delete m_A;
      m_A = m_diffEq.createMatrix();
      delete m_solver;
      m_solver = SolverUtils::instance()->createLinearSolver(m_A);
    }
  }

  void setZero()
  {
    x0.setZero();
    m_rhs.setZero();
    m_lhs.setZero();
    m_increment.setZero();
    m_A->setZero();
  }

  Scalar computeResidual()
  {
    m_rhs.setZero();
    m_diffEq.evaluatePDot(m_rhs);

    m_lhs.setZero();
    for (size_t i = 0; i < m_fixed.size(); ++i) {
      int idx = m_fixed[i];
      m_rhs[idx] = m_desired[i];
      m_lhs[idx] = m_diffEq.getX(idx);
    }

    m_infnorm = (m_lhs - m_rhs).lpNorm<Eigen::Infinity>();

    return (m_lhs - m_rhs).norm();
  }

  bool isConverged()
  {
    m_residual = computeResidual();
//     std::cout << "atol " << m_residual << std::endl
//               << "infnorm " << m_infnorm << std::endl
//               << "rtol " << m_residual / m_initial_residual << std::endl
//               << "stol " << m_increment.norm() << std::endl;
    if ( m_residual < m_atol ) {
      //std::cout << "converged atol" << std::endl;
      return true;
    }
    if ( m_infnorm < m_inftol ) {
      //std::cout << "converged inftol" << std::endl;
      return true;
    }
    if ( m_residual / m_initial_residual < m_rtol ) {
      //std::cout << "converged rtol" << std::endl;
      return true;
    }
    if ( m_increment.norm() < m_stol ) {
      //std::cout << "converged stol" << std::endl;
      return true;
    }
    return false;
  }

protected:

  ODE& m_diffEq;

  VecXd x0;
  VecXd m_rhs;
  VecXd m_lhs;
  VecXd m_increment;

  IntArray m_fixed;
  std::vector<Scalar> m_desired;

  Scalar m_initial_residual;
  Scalar m_residual;
  int m_curit;

  MatrixBase* m_A;
  LinearSolverBase* m_solver;
};

} // namespace BASim

#endif // STATICSSOLVER_HH

/**
 * \file DiffEqSolver.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef DIFFEQSOLVER_HH
#define DIFFEQSOLVER_HH

#include "BASim/src/Core/Definitions.hh"

namespace BASim {

/** Base class for all time steppers. */
class DiffEqSolver
{
public:

  virtual ~DiffEqSolver() {}

  Scalar getTime() const { return m_time; }

  void setTime(Scalar time) { m_time = time; }

  Scalar getTimeStep() const { return m_dt; }

  void setTimeStep(Scalar dt) { m_dt = dt;} 

  virtual void execute() = 0;

  virtual std::string getName() const = 0;

  int getMaxIterations() const { return m_maxit; }

  void setMaxIterations(int iterations)
  {
    m_maxit = std::max(1, iterations);
  }

  Scalar get_stol() const { return m_stol; }
  void set_stol(Scalar s) { m_stol = s; }

  Scalar get_atol() const { return m_atol; }
  void set_atol(Scalar a) { m_atol = a; }

  Scalar get_rtol() const { return m_rtol; }
  void set_rtol(Scalar r) { m_rtol = r; }

  Scalar get_inftol() const { return m_inftol; }
  void set_inftol(Scalar i) { m_inftol = i; }

protected:

  DiffEqSolver(Scalar time = 0, Scalar dt = 0.1)
    : m_time(time)
    , m_dt(dt)
    , m_maxit(50)
    , m_stol(1e-8)
    , m_atol(1e-8)
    , m_rtol(1e-8)
    , m_inftol(1e-8)
    , m_infnorm(0)
  {}

  Scalar m_time; ///< the current time
  Scalar m_dt; ///< size of the time step

  int m_maxit; ///< maximum number of iterations (only for implicit methods)

  Scalar m_stol; ///< minimum change in norm of solution between iterations
  Scalar m_atol; ///< tolerance for absolute norm of solution
  Scalar m_rtol; ///< tolerance for norm of solution relative to relative to initial guess
  Scalar m_inftol; ///< tolerance for the infinity norm of the residual
  Scalar m_infnorm;
};

} // namespace BASim

#endif // DIFFEQSOLVER_HH

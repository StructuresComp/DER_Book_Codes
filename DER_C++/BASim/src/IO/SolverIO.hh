/**
 * \file SolverIO.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 11/17/2009
 */

#ifndef SOLVERIO_HH
#define SOLVERIO_HH

#include "BASim/src/Math/LinearSolverBase.hh"
#include "BASim/src/Math/SolverUtils.hh"

namespace BASim {

#define CHECK_RETURN(check) { if (check != 0) { ret = -1; break; } }

inline int parseSolver(Tokenizer& tokenizer)
{
  std::string solver;
  int ret = tokenizer.read(solver);
  if (ret == -1) {
    std::cerr << "Error parsing solver type" << std::endl;
    return -1;
  }

  if (solver == "auto") {
    SolverUtils::instance()->setSolverType(SolverUtils::AUTO_SOLVER);

  } else if (solver == "conjugate-gradient") {
    SolverUtils::instance()->setSolverType(SolverUtils::CONJUGATE_GRADIENT);

#ifdef HAVE_PARDISO
  } else if (solver == "pardiso") {
    SolverUtils::instance()->setSolverType(SolverUtils::PARDISO_SOLVER);
#endif // HAVE_PARDISO

#ifdef HAVE_PETSC
  } else if (solver == "petsc") {
    SolverUtils::instance()->setSolverType(SolverUtils::PETSC_SOLVER);
#endif // HAVE_PETSC

#ifdef HAVE_MKL
  } else if (solver == "mkl") {
    SolverUtils::instance()->setSolverType(SolverUtils::MKL_LINEAR_SOLVER);
#endif // HAVE_MKL

#ifdef HAVE_LAPACK
  } else if (solver == "lapack") {
    SolverUtils::instance()->setSolverType(SolverUtils::MKL_LINEAR_SOLVER);
#endif // HAVE_LAPACK

  } else {
    std::cerr << "Unknown solver type " << solver << std::endl;
    return -1;
  }

  return 0;
}

inline int parseMatrix(Tokenizer& tokenizer)
{
  std::string matrix;
  int ret = tokenizer.read(matrix);
  if (ret == -1) {
    std::cerr << "Error parsing matrix type" << std::endl;
    return -1;
  }

  if (matrix == "auto") {
    SolverUtils::instance()->setMatrixType(SolverUtils::AUTO_MATRIX);

  } else if (matrix == "band-matrix") {
    SolverUtils::instance()->setMatrixType(SolverUtils::BAND_MATRIX);

#ifdef HAVE_PARDISO
  } else if (matrix == "pardiso") {
    SolverUtils::instance()->setMatrixType(SolverUtils::PARDISO_MATRIX);
#endif // HAVE_PARDISO

#ifdef HAVE_PETSC
  } else if (matrix == "petsc") {
    SolverUtils::instance()->setMatrixType(SolverUtils::PETSC_MATRIX);
#endif // HAVE_PETSC

  } else {
    std::cerr << "Unknown matrix type " << matrix << std::endl;
    return -1;
  }

  return 0;
}

inline int readSolverFile(const std::string& file)
{
  int ret = 0;

  std::ifstream in_file;
  in_file.open(file.c_str());
  if (!in_file.is_open()) {
    std::cerr << "Failed to open file " << file << std::endl;
    return -1;
  }

  Tokenizer tokenizer(in_file);
  std::string token;

  while (tokenizer.read(token) == 0) {
    if (token == "solver-type") {
      int check = parseSolver(tokenizer);
      CHECK_RETURN(check);
    } else if (token == "matrix-type") {
      int check = parseMatrix(tokenizer);
      CHECK_RETURN(check);
    } else {
      std::cerr << "Unknown token in file " << file << ": " << token
                << std::endl;
      ret = -1;
      break;
    }
  }

  in_file.close();
  return ret;
}

}

#endif // SOLVERIO_HH

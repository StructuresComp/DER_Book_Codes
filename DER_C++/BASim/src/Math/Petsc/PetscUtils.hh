/**
 * \file PetscUtils.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/08/2009
 */

#ifndef PETSCUTILS_HH
#define PETSCUTILS_HH

#include <petscmat.h>

namespace BASim {

namespace PetscUtils {

inline void initializePetsc(int* argc, char*** argv)
{
  PetscInitialize(argc, argv, PETSC_NULL, PETSC_NULL);
  PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
}

inline void finalizePetsc()
{
  PetscFinalize();
}

template <class Vector>
inline int copyToPetscVector(Vec petscVec, const Vector& otherVec)
{
  PetscScalar* petscDat;
  int ierr = VecGetArray(petscVec, &petscDat);
  CHKERRQ(ierr);
  for (int i = 0; i < otherVec.size(); ++i) petscDat[i] = otherVec[i];
  ierr = VecRestoreArray(petscVec, &petscDat);
  CHKERRQ(ierr);
  return 0;
}

inline int createPetscVector(Vec& petscVec, PetscInt n)
{
  if (petscVec != NULL) return 0;
  int ierr = VecCreateSeq(PETSC_COMM_SELF, n, &petscVec);
  CHKERRQ(ierr);
  return 0;
}

template <class Vector>
inline int copyFromPetscVector(Vector& otherVec, const Vec& petscVec)
{
  PetscScalar* petscDat;
  int ierr = VecGetArray(const_cast<Vec&>(petscVec), &petscDat);
  CHKERRQ(ierr);
  for (int i = 0; i < otherVec.size(); ++i) otherVec[i] = petscDat[i];
  ierr = VecRestoreArray(petscVec, &petscDat);
  CHKERRQ(ierr);
  return 0;
}

template <class Vector>
inline int addToPetscVector(Vec& petscVec, Scalar s, const Vector& otherVec)
{
  PetscScalar* petscDat;
  int ierr = VecGetArray(petscVec, &petscDat);
  CHKERRQ(ierr);
  for (int i = 0; i < otherVec.size(); ++i) petscDat[i] += s * otherVec[i];
  ierr = VecRestoreArray(petscVec, &petscDat);
  CHKERRQ(ierr);
  return 0;
}

template <class Vector>
inline int addFromPetscVector(Vector& otherVec, Scalar s, const Vec& petscVec)
{
  PetscScalar* petscDat;
  int ierr = VecGetArray(petscVec, &petscDat);
  CHKERRQ(ierr);
  for (int i = 0; i < otherVec.size(); ++i) otherVec[i] += s * petscDat[i];
  ierr = VecRestoreArray(petscVec, &petscDat);
  CHKERRQ(ierr);
  return 0;
}

} // namespace PetscUtils

} // namespace BASim

#endif // PETSCUTILS_HH

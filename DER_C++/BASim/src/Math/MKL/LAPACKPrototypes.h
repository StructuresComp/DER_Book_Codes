/**
 * \file LAPACKPrototypes.h
 *
 * \author miklos@cs.columbia.edu
 * \date 03/11/2010
 */

#ifndef LAPACKPROTOTYPES_H
#define LAPACKPROTOTYPES_H

#ifdef HAVE_MKL
#include <mkl_lapack.h>
#else // HAVE_MKL
extern "C" void dgbsv_( int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info );
#endif // HAVE_MKL

#endif // LAPACKPROTOTYPES_H

/**
 * \file EigenIncludes.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/15/2009
 */

#ifndef EIGENINCLUDES_HH
#define EIGENINCLUDES_HH

//#define EIGEN_DONT_ALIGN

#define EIGEN_VECTOR_IO Eigen::IOFormat(8, Eigen::Raw, ",", ",", "", "", "{", "}")
#define EIGEN_MATRIX_IO Eigen::IOFormat(8, Eigen::Raw, ",", ",", "{", "}", "{", "}")
#define EIGEN_SPACES_ONLY_IO Eigen::IOFormat(8, Eigen::Raw, " ", " ", "", "", "", "")
#define EIGEN_DEFAULT_IO_FORMAT EIGEN_MATRIX_IO

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO 1

// Khalid Jawed
// We were having trouble locating Eigen. So, Eigen is now included in
// "src" folder and linked by hard code (below)
#include <BASim/src/Eigen2/Eigen/LU>
#include <BASim/src/Eigen2/Eigen/Cholesky>
#include <BASim/src/Eigen2/Eigen/QR>
#include <BASim/src/Eigen2/Eigen/SVD>
#include <BASim/src/Eigen2/Eigen/Dense>
#include <BASim/src/Eigen2/Eigen/Core>
#include <BASim/src/Eigen2/Eigen/Geometry>
#include <BASim/src/Eigen2/Eigen/StdVector>

#endif // EIGENINCLUDES_HH

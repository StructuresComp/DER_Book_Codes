cmake_minimum_required (VERSION 2.6)

project (BASim)
set (BASIM_MAJOR_VERSION "0")
set (BASIM_MINOR_VERSION "1")
set (BASIM_BUGFIX_VERSION "0")
set (BASIM_VERSION_NUMBER "${BASIM_MAJOR_VERSION}.${BASIM_MINOR_VERSION}.${BASIM_BUGFIX_VERSION}")
set (BASIM_VERSION "${BASIM_VERSION_NUMBER}")
add_definitions (-DBASIM_VERSION="${BASIM_VERSION}")

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Debug CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# add directory with macros
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# add cmake file
include (BASimCommon)

# set options
option (BUILD_APPS "Build applications" ON)
option (TEST_ROD_STRETCHING_FORCE "Tests RodStretchingForce" OFF)
option (TEST_ROD_TWISTING_FORCE "Tests RodTwistingForce" OFF)
option (TEST_ROD_BENDING_FORCE "Tests RodBendingForce" OFF)
option (TIMING_ON "Enables timing macros" OFF)
if (TIMING_ON)
  add_definitions (-DTIMING_ON)
endif (TIMING_ON)

# Add some extra warnings, sse support
if (CMAKE_COMPILER_IS_GNUCXX)
  # LINUX
  if (CMAKE_SYSTEM_NAME MATCHES Linux)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor -ansi -Wcast-align -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wformat-security -fno-common -Wno-unused -Wno-empty-body")
    if (CMAKE_BUILD_TYPE MATCHES Release)
      add_definitions (-DNDEBUG)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mfpmath=sse -msse2 -fno-check-new -fstrict-aliasing") # -fno-exceptions
    endif (CMAKE_BUILD_TYPE MATCHES Release)
  endif (CMAKE_SYSTEM_NAME MATCHES Linux)
  # OS X
  if (CMAKE_SYSTEM_NAME MATCHES Darwin)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor -ansi -Wcast-align -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wformat-security -fno-common -Wno-unused -Wno-empty-body")
    if (CMAKE_BUILD_TYPE MATCHES Release)
      add_definitions (-DNDEBUG)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse -msse2 -fno-check-new -fstrict-aliasing") # -fno-exceptions
    endif (CMAKE_BUILD_TYPE MATCHES Release)
  endif (CMAKE_SYSTEM_NAME MATCHES Darwin)
endif (CMAKE_COMPILER_IS_GNUCXX)

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions (-DDEBUG)
endif (CMAKE_BUILD_TYPE MATCHES Debug)

# Eigen library is required
find_package (Eigen REQUIRED)
include_directories (${EIGEN_INCLUDE_DIR})

option (USE_PNG "Builds in support for exporting to png" ON)
if (USE_PNG)
  find_package (PNG)
  if (PNG_FOUND)
    add_definitions (-DHAVE_PNG)
    add_definitions (${PNG_DEFINITIONS})
    include_directories (${PNG_INCLUDE_DIR})
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${PNG_LIBRARIES})
  endif (PNG_FOUND)
endif (USE_PNG)

option (USE_PETSc "Builds in support for PETSc" OFF)
if (USE_PETSc)
  find_package (PETSc)
  if (PETSC_FOUND)
    add_definitions (-DHAVE_PETSC)
    include_directories (${PETSC_INCLUDE_PATH})
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${PETSC_LIBRARIES})
  endif (PETSC_FOUND)
endif (USE_PETSc)

#find_package (FLENS REQUIRED)
#include_directories (${FLENS_INCLUDE_DIR})

option (USE_MKL "Builds in support for Intel's MKL" ON)
if (USE_MKL)
  find_package (MKL)
  if (MKL_FOUND)
    add_definitions (-DHAVE_LAPACK -DHAVE_MKL)
    include_directories (${MKL_INCLUDES})
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${MKL_LIBRARIES})
  else (MKL_FOUND)
    message ("MKL not found. Try using generic LAPACK instead by setting USE_LAPACK to ON")
  endif (MKL_FOUND)
endif (USE_MKL)

option (USE_LAPACK "Builds in support for generic LAPACK if MKL is not available" OFF)
if (USE_LAPACK AND NOT MKL_FOUND)
  find_package (LAPACK)
  if (LAPACK_FOUND)
    add_definitions (-DHAVE_LAPACK)
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${LAPACK_LIBRARIES})
  endif (LAPACK_FOUND)
endif (USE_LAPACK AND NOT MKL_FOUND)

option (USE_PARDISO "Builds in support for Pardiso" OFF)
if (USE_PARDISO)
  find_package (PARDISO)
  if (PARDISO_FOUND)
    add_definitions (-DHAVE_PARDISO)
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${PARDISO_LIBRARIES})
  endif (PARDISO_FOUND)
endif (USE_PARDISO)

option (USE_OPENMP "Builds in support for multi-threading using OpenMP" OFF)
if (USE_OPENMP)
  find_package (OpenMP)
  if (OPENMP_FOUND)
    add_definitions (-DHAVE_OPENMP)
    #set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif (OPENMP_FOUND)
endif (USE_OPENMP)


include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory (BASim)

if (BUILD_APPS)
  add_subdirectory (Apps)
endif (BUILD_APPS)

add_subdirectory (doc)

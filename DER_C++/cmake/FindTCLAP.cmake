# - Find TCLAP
# Find the TCLAP headers
#
# TCLAP_INCLUDE_DIR - where to find the TCLAP headers
# TCLAP_FOUND       - True if TCLAP is found

if (TCLAP_INCLUDE_DIR)
  # already in cache, be silent
  set (TCLAP_FIND_QUIETLY TRUE)
endif (TCLAP_INCLUDE_DIR)

# find the headers
find_path (TCLAP_INCLUDE_PATH
  tclap/CmdLine.h
  HINTS TCLAP_INC_DIR ENV TCLAP_INC_DIR
  DOC "TCLAP include directory (http://tclap.sourceforge.net/)"
  )

# handle the QUIETLY and REQUIRED arguments and set TCLAP_FOUND to
# TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (TCLAP "TCLAP (http://tclap.sourceforge.net/) could not be found. Set TCLAP_INC_DIR to point to the headers." TCLAP_INCLUDE_PATH)

if (TCLAP_FOUND)
  set (TCLAP_INCLUDE_DIR ${TCLAP_INCLUDE_PATH})
endif (TCLAP_FOUND)

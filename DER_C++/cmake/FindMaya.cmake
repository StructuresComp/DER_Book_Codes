# - Find Maya
# Find the Maya headers and libraries
#
# MAYA_INCLUDE_DIR - include path for headers
# MAYA_LIBRARIES   - libraries to include when linking
# MAYA_FOUND       - True if MAYA is found

if (MAYA_INCLUDE_DIR AND MAYA_LIBRARIES)
  # already in cache, be silent
  set (MAYA_FIND_QUIETLY TRUE)
endif (MAYA_INCLUDE_DIR AND MAYA_LIBRARIES)


##############################################################
##
## Search for the header location. There are hundreds of
## headers so just search for one we know must be there and 
## assume it is a valid install.
##
##############################################################
find_path (MAYA_INCLUDE_PATH
  /include/maya/MPxNode.h
  HINTS ENV MAYA_LOCATION
  )

##############################################################
##
## Search for all the libraries to make sure they can be found
##
##############################################################
set (MAYA_LIBRARIES_FOUND)
set (MAYA_LIBRARIES_MISSING)
set (MAYA_LIBS OpenMaya OpenMayaUI OpenMayaAnim OpenMayaFX OpenMayaRender)
foreach (MAYALIB ${MAYA_LIBS})
    set (MAYA_SEARCH_LIB "MAYA_SEARCH_LIB-NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
    find_library (MAYA_SEARCH_LIB ${MAYALIB} PATHS $ENV{MAYA_LOCATION}/lib)
    if (MAYA_SEARCH_LIB)
        list (APPEND MAYA_LIBRARIES_FOUND ${MAYA_SEARCH_LIB})
    else (MAYA_SEARCH_LIB)
        list (APPEND MAYA_LIBRARIES_MISSING ${MAYA_SEARCH_LIB})
        message (SEND_ERROR "Unable to find Maya library ${MAYALIB}")
    endif (MAYA_SEARCH_LIB)
endforeach (MAYALIB)

set (MAYA_LIBRARY ${MAYA_LIBRARIES_FOUND} CACHE STRING "Maya libraries" FORCE)

# handle the QUIETLY and REQUIRED arguments and set
# MAYA_FOUND to TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
#find_package_handle_standard_args (MAYA "MAYA not found. Set MAYA_LOCATION " MAYA_INCLUDE_PATH MAYA_LIBRARIES)
find_package_handle_standard_args (MAYA "MAYA not found. Set MAYA_LOCATION " MAYA_INCLUDE_PATH MAYA_LIBRARY)

if (MAYA_FOUND)
  set (MAYA_INCLUDE_DIR ${MAYA_INCLUDE_PATH}/include)
  set (MAYA_LIBRARIES ${MAYA_LIBRARY})
  set (MAYA_CXX_FLAGS "-Dlinux -DLINUX -DREQUIRE_IOSTREAM -D_BOOL -fPIC -g")
else (MAYA_FOUND)
  set (MAYA_INCLUDE_DIR)
  set (MAYA_LIBRARIES)
endif (MAYA_FOUND)

mark_as_advanced (MAYA_LIBRARIES MAYA_SEARCH_LIB MAYA_LIBRARY)

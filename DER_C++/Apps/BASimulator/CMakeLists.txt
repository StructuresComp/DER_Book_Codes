append_files (Headers "hh" . Problems)
append_files (Headers "inl" . Problems)
append_files (Templates "tcc" . Problems)
append_files (Sources "cc" . Problems)

find_package (TCLAP REQUIRED)
include_directories (${TCLAP_INCLUDE_DIR})

#find_package (wxWidgets REQUIRED base core gl)
#include (${wxWidgets_USE_FILE})

#set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} -lnoise)
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES})

add_executable (BASimulator ${Headers} ${Templates} ${Sources})
#target_link_libraries (BASimulator BASim -lnoise)
# target_link_libraries (BASimulator BASim -lfftw3)
target_link_libraries (BASimulator BASim)

# Creates a link to the assets directory in the binary output directory so the binary can access obj files and such
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/assets ${CMAKE_CURRENT_BINARY_DIR}/assets )

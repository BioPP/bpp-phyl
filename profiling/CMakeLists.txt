# Authors:
#   Francois Gindraud (2017)
# Created: 23/06/2017

# Script that compiles a local version of bpp-phyl, which separates NewPhyl from the rest.

set (PROF_BPP_SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set (PROF_BPP_BIN_DIR ${CMAKE_BINARY_DIR}/src)

# Compile a lib with only Bpp/Phyl/ stuff ; standard compile flags
file (GLOB_RECURSE PROF_PHYL_CPP_FILES ${PROF_BPP_SRC_DIR}/Bpp/Phyl/*.cpp)
add_library (prof-lib-legacy STATIC EXCLUDE_FROM_ALL ${PROF_PHYL_CPP_FILES})
target_include_directories (prof-lib-legacy PUBLIC
  $<BUILD_INTERFACE:${PROF_BPP_SRC_DIR}>
  $<BUILD_INTERFACE:${PROF_BPP_BIN_DIR}>
  $<INSTALL_INTERFACE:include>
  )
target_link_libraries (prof-lib-legacy ${BPP_LIBS_STATIC} Eigen3::Eigen)
target_compile_options (prof-lib-legacy
  PRIVATE
  -O2
  -g
  -fno-omit-frame-pointer
  )

# Compile a program with custom debug flags
file (GLOB_RECURSE PROF_NEWPHYL_CPP_FILES ${PROF_BPP_SRC_DIR}/Bpp/NewPhyl/*.cpp)
set (PROF_EXEC_TARGET guinea-pig)

add_executable (${PROF_EXEC_TARGET} EXCLUDE_FROM_ALL
  ${PROF_NEWPHYL_CPP_FILES}
  ${CMAKE_SOURCE_DIR}/test/new_phyl_dataflow.cpp
  )
target_link_libraries (${PROF_EXEC_TARGET} ${BPP_LIBS_STATIC} prof-lib-legacy)
target_compile_options (${PROF_EXEC_TARGET}
  PRIVATE
  -O2
  -g
  -fno-omit-frame-pointer
  #-fno-inline-functions
  #gcc:  -fno-inline-functions-called-once
  #-fno-optimize-sibling-calls
  )

# Convenience targets
add_custom_target (
  prof-perf
  COMMAND perf record -g -- ${CMAKE_CURRENT_BINARY_DIR}/${PROF_EXEC_TARGET}
  DEPENDS ${PROF_EXEC_TARGET}
  COMMENT "Perf record on ${PROF_EXEC_TARGET}"
  )


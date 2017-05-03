# CMake script for Bio++ PhylLib
# Authors:
#   Sylvain Gaillard
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 20/08/2009

cmake_minimum_required (VERSION 2.8.12)
project (bpp-phyl CXX)

# Compile options
set (private-compile-options -std=c++11 -Wall -Weffc++ -Wshadow -Wconversion)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(NOT NO_DEP_CHECK)
  SET(NO_DEP_CHECK FALSE CACHE BOOL
      "Disable dependencies check for building distribution only."
      FORCE)
ENDIF(NOT NO_DEP_CHECK)

IF(NO_DEP_CHECK)
  MESSAGE("-- Dependencies checking disabled. Only distribution can be built.")
ELSE(NO_DEP_CHECK)

# Libtool-like version number
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current:  The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age:      The difference between the newest and oldest interfaces that this
#           library implements.
# In other words, the library implements all the interface numbers in the
# range from number current - age to current.
SET(${PROJECT_NAME}_VERSION_CURRENT "10")
SET(${PROJECT_NAME}_VERSION_REVISION "3")
SET(${PROJECT_NAME}_VERSION_AGE "1")

# Effective version number computation
MATH(EXPR ${PROJECT_NAME}_VERSION_MAJOR "${${PROJECT_NAME}_VERSION_CURRENT} - ${${PROJECT_NAME}_VERSION_AGE}")
SET(${PROJECT_NAME}_VERSION_MINOR ${${PROJECT_NAME}_VERSION_AGE})
SET(${PROJECT_NAME}_VERSION_PATCH ${${PROJECT_NAME}_VERSION_REVISION})
SET(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

set (PROJECT_VERSION ${${PROJECT_NAME}_VERSION})

# Find dependencies (add install directory to search)
if (CMAKE_INSTALL_PREFIX)
  set (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif (CMAKE_INSTALL_PREFIX)

find_package (bpp-core 2.0.3 REQUIRED)
find_package (bpp-seq 9.1.3 REQUIRED)

# CMake package
set (cmake-package-location lib/cmake/${PROJECT_NAME})
include (CMakePackageConfigHelpers)
configure_package_config_file (
  package.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/package.cmake
  INSTALL_DESTINATION ${cmake-package-location}
  )
write_basic_package_version_file (
  ${CMAKE_CURRENT_BINARY_DIR}/package-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  )
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/package.cmake DESTINATION ${cmake-package-location}
  RENAME ${PROJECT_NAME}-config.cmake)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/package-version.cmake DESTINATION ${cmake-package-location}
  RENAME ${PROJECT_NAME}-config-version.cmake)

# Check demangling support
include (CheckCXXSymbolExists)
CHECK_CXX_SYMBOL_EXISTS("abi::__cxa_demangle" "cxxabi.h" BPP_HAVE_DEMANGLING)

# Use eigen (TODO, add eigen files to repository directly)
find_package (Eigen3)

# Define the libraries
add_subdirectory (src)

# Doxygen
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND)
  ADD_CUSTOM_TARGET (apidoc cp Doxyfile ${CMAKE_BINARY_DIR}/Doxyfile-build
    COMMAND echo "OUTPUT_DIRECTORY=${CMAKE_BINARY_DIR}" >> ${CMAKE_BINARY_DIR}/Doxyfile-build
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile-build
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  ADD_CUSTOM_TARGET (apidoc-stable cp Doxyfile ${CMAKE_BINARY_DIR}/Doxyfile-stable
    COMMAND echo "OUTPUT_DIRECTORY=${CMAKE_BINARY_DIR}" >> ${CMAKE_BINARY_DIR}/Doxyfile-stable
    COMMAND echo "HTML_HEADER=header.html" >> ${CMAKE_BINARY_DIR}/Doxyfile-stable
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile-stable
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
ENDIF (DOXYGEN_FOUND)

ENDIF(NO_DEP_CHECK)

# Packager
SET(CPACK_PACKAGE_NAME "libbpp-phyl")
SET(CPACK_PACKAGE_VENDOR "Bio++ Development Team")
SET(CPACK_PACKAGE_VERSION "2.2.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "2")
SET(CPACK_PACKAGE_VERSION_MINOR "2")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Bio++ Phylogenetics library")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES
 "CMakeFiles"
 "Makefile"
 "_CPack_Packages"
 "CMakeCache.txt"
 ".*\\\\.cmake"
 ".*\\\\.git"
 ".*\\\\.gz"
 ".*\\\\.deb"
 ".*\\\\.rpm"
 ".*\\\\.dmg"
 ".*\\\\.sh"
 ".*\\\\..*\\\\.swp"
 "src/\\\\..*"
 "src/libbpp*"
 "test/tmp_*"
 "tmp_*"
 "html"
 "Phyl.tag"
 "Testing"
 "build-stamp"
 "install_manifest.txt"
 "DartConfiguration.tcl"
 ".clang_format"
 ".clang_complete"
 ${CPACK_SOURCE_IGNORE_FILES}
)
IF (MACOS)
  SET(CPACK_GENERATOR "Bundle")
ENDIF()

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "lib${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
# 'clean' is not (yet) a first class target. However, we need to clean the directories before building the sources:
IF("${CMAKE_GENERATOR}" MATCHES "Make")
  ADD_CUSTOM_TARGET(make_clean
  COMMAND ${CMAKE_MAKE_PROGRAM} clean
  WORKING_DIRECTORY ${CMAKE_CURRENT_DIR}
  )
  ADD_DEPENDENCIES(dist make_clean)
ENDIF()

IF(NOT NO_DEP_CHECK)
IF (UNIX)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist)

ENDIF()

SET(CTEST_UPDATE_TYPE git)
SET(UPDATE_COMMAND "git")
SET(UPDATE_OPTIONS "")

enable_testing ()
include (CTest)
if (BUILD_TESTING)
  add_subdirectory (test)
endif ()

ENDIF(NOT NO_DEP_CHECK)

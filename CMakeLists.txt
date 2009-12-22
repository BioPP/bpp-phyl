# CMake script for Bio++ PhylLib
# Author: Sylvain Gaillard and Julien Dutheil
# Created: 20/08/2009

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(bppphyl CXX)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_CXX_FLAGS "-Wall -Wshadow")
IF(NOT NO_VIRTUAL_COV)
  SET(NO_VIRTUAL_COV FALSE CACHE BOOL
      "Disable covariant return type with virtual inheritance, for compilers that do not support it."
      FORCE)
ENDIF(NOT NO_VIRTUAL_COV)

IF(NO_VIRTUAL_COV)
  MESSAGE("-- Covariant return with virtual inheritance disabled.")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_VIRTUAL_COV=1")
ENDIF(NO_VIRTUAL_COV)

# Libtool-like version number
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current:  The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age:      The difference between the newest and oldest interfaces that this
#           library implements.
# In other words, the library implements all the interface numbers in the
# range from number current - age to current.
SET(BPPPHYL_VERSION_CURRENT "8")
SET(BPPPHYL_VERSION_REVISION "0")
SET(BPPPHYL_VERSION_AGE "0")

# Effective version number computation
MATH(EXPR BPPPHYL_VERSION_MAJOR "${BPPPHYL_VERSION_CURRENT} - ${BPPPHYL_VERSION_AGE}")
SET(BPPPHYL_VERSION_MINOR ${BPPPHYL_VERSION_AGE})
SET(BPPPHYL_VERSION_PATCH ${BPPPHYL_VERSION_REVISION})
SET(BPPPHYL_VERSION "${BPPPHYL_VERSION_MAJOR}.${BPPPHYL_VERSION_MINOR}.${BPPPHYL_VERSION_PATCH}")

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

#here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  SET(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
  FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})

  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} found here:")
    MESSAGE("   includes : ${${lib_name}_INCLUDE_DIR}")
    MESSAGE("   libraries: ${${lib_name}_LIBRARY}")
  ELSE(${lib_name}_LIBRARY)
    MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
  ENDIF(${lib_name}_LIBRARY)
  
  #add the dependency:
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(IMPROVED_FIND_LIBRARY)

#Find the Bio++ libraries:
IMPROVED_FIND_LIBRARY(LIBS bpp-seq Seq/Alphabet.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-numcalc NumCalc/NumTools.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-utils Utils/Clonable.h)

# Subdirectories
ADD_SUBDIRECTORY(Phyl)

# Doxygen
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND)
  ADD_CUSTOM_TARGET (apidoc ${DOXYGEN_EXECUTABLE} Doxyfile)
ENDIF (DOXYGEN_FOUND)

# Packager
SET(CPACK_PACKAGE_NAME "libbpp-phyl")
SET(CPACK_PACKAGE_VENDOR "Bio++ Development Team")
SET(CPACK_PACKAGE_VERSION "1.8.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "8")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Bio++ Phylogenetics library")
SET(CPACK_RESOURCE_FILE_LICENSE "COPYING")
SET(CPACK_RESOURCE_FILE_AUTHORS "AUTHORS")
SET(CPACK_GENERATOR "DEB;RPM")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SYSTEM_NAME x86_64)
SET(CPACK_TOPLEVEL_TAG x86_64)
SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64)
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6, libbpp-utils (= 1.4.0), libbpp-numcalc (= 1.7.0), libbpp-seq (= 1.6.0)")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Julien Dutheil <jdutheil@birc.au.dk>")
SET(CPACK_SOURCE_IGNORE_FILES
 ".*CMakeFiles"
 ".*Makefile"
 ".*_CPack_Packages"
 ".*CMakeCache.txt"
 "debian/"
 "Phyl.spec"
 ".*\\\\.cmake"
 ".*CVS/"
 ".*.gz"
 ".*.deb"
 ".*.rpm"
 ".*~"
 "Phyl/\\\\..*"
 "html"
 "Phyl.tag"
 "autogen.sh"
 "autoclean.sh"
 ${CPACK_SOURCE_IGNORE_FILES}
)
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)


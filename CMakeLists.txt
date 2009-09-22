# CMake script for Bio++ Utils
# Author: Sylvain Gaillard
# Created: 20/08/2009

# Global parameters
cmake_minimum_required(VERSION 2.6)
project(bppphyl CXX)
if(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT DEFINED CMAKE_BUILD_TYPE)
set(CMAKE_CXX_FLAGS "-Wall")

# Libtool-like version number
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current:  The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age:      The difference between the newest and oldest interfaces that this
#           library implements.
# In other words, the library implements all the interface numbers in the
# range from number current - age to current.
set(BPPPHYL_VERSION_CURRENT "8")
set(BPPPHYL_VERSION_REVISION "0")
set(BPPPHYL_VERSION_AGE "0")

# Effective version number computation
math(EXPR BPPPHYL_VERSION_MAJOR "${BPPPHYL_VERSION_CURRENT} - ${BPPPHYL_VERSION_AGE}")
set(BPPPHYL_VERSION_MINOR ${BPPPHYL_VERSION_AGE})
set(BPPPHYL_VERSION_PATCH ${BPPPHYL_VERSION_REVISION})
set(BPPPHYL_VERSION "${BPPPHYL_VERSION_MAJOR}.${BPPPHYL_VERSION_MINOR}.${BPPPHYL_VERSION_PATCH}")

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
if(CMAKE_INSTALL_PREFIX)
  set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif(CMAKE_INSTALL_PREFIX)

# Find the bpp-utils library
find_path(BPPUTILS_INCLUDE_DIR Utils/Clonable.h)
set(BPPUTILS_NAMES bpp-utils bpp-utilslib bpp-utilsdll)
find_library(BPPUTILS_LIBRARY NAMES ${BPPUTILS_NAMES})
set( BPPUTILS_LIBRARIES ${BPPUTILS_LIBRARY} )
if(BPPUTILS_LIBRARY)
  message("-- bpp-utils found here:")
  message("   includes : ${BPPUTILS_INCLUDE_DIR}")
  message("   libraries: ${BPPUTILS_LIBRARIES}")
else(BPPUTILS_LIBRARY)
  message(FATAL_ERROR "bpp-utils required but not found")
endif(BPPUTILS_LIBRARY)

# Include the bpp-utils dependency
include_directories(${BPPUTILS_INCLUDE_DIR})
set(LIBS ${LIBS} ${BPPUTILS_LIBRARIES})

# Find the bpp-numcalc library
find_path(BPPNUMCALC_INCLUDE_DIR NumCalc/NumTools.h)
set(BPPNUMCALC_NAMES bpp-numcalc bpp-numcalclib bpp-numcalcdll)
find_library(BPPNUMCALC_LIBRARY NAMES ${BPPNUMCALC_NAMES})
set( BPPNUMCALC_LIBRARIES ${BPPNUMCALC_LIBRARY} )
if(BPPNUMCALC_LIBRARY)
  message("-- bpp-numcalc found here:")
  message("   includes : ${BPPNUMCALC_INCLUDE_DIR}")
  message("   libraries: ${BPPNUMCALC_LIBRARIES}")
else(BPPNUMCALC_LIBRARY)
  message(FATAL_ERROR "bpp-numcalc required but not found")
endif(BPPNUMCALC_LIBRARY)

# Include the bpp-numcalc dependency
include_directories(${BPPNUMCALC_INCLUDE_DIR})
set(LIBS ${LIBS} ${BPPNUMCALC_LIBRARIES})

# Find the bpp-seq library
find_path(BPPSEQ_INCLUDE_DIR Seq/Alphabet.h)
set(BPPSEQ_NAMES bpp-seq bpp-seqlib bpp-seqdll)
find_library(BPPSEQ_LIBRARY NAMES ${BPPSEQ_NAMES})
set( BPPSEQ_LIBRARIES ${BPPSEQ_LIBRARY} )
if(BPPSEQ_LIBRARY)
  message("-- bpp-seq found here:")
  message("   includes : ${BPPSEQ_INCLUDE_DIR}")
  message("   libraries: ${BPPSEQ_LIBRARIES}")
else(BPPSEQ_LIBRARY)
  message(FATAL_ERROR "bpp-seq required but not found")
endif(BPPSEQ_LIBRARY)

# Include the bpp-seq dependency
include_directories(${BPPSEQ_INCLUDE_DIR})
set(LIBS ${LIBS} ${BPPSEQ_LIBRARIES})

# Subdirectories
add_subdirectory(Phyl)

# Doxygen
find_package(Doxygen)
if (DOXYGEN_FOUND)
  add_custom_target (apidoc ${DOXYGEN_EXECUTABLE} Doxyfile)
endif (DOXYGEN_FOUND)

# Packager
set(CPACK_PACKAGE_NAME "libbpp-phyl")
set(CPACK_PACKAGE_VENDOR "Bio++ Development Team")
set(CPACK_PACKAGE_VERSION "1.8.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "8")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Bio++ Phylogenetics library")
set(CPACK_RESOURCE_FILE_LICENSE "COPYING")
set(CPACK_RESOURCE_FILE_AUTHORS "AUTHORS")
set(CPACK_GENERATOR "DEB;RPM")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SYSTEM_NAME x86_64)
set(CPACK_TOPLEVEL_TAG x86_64)
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE amd64)
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6, libbpp-utils (= 1.4.0), libbpp-numcalc (= 1.7.0), libbpp-seq (= 1.6.0)")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Julien Dutheil <jdutheil@birc.au.dk>")
set(CPACK_SOURCE_IGNORE_FILES
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
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
include(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)


# CMake script for bpp-phyl unit tests
# Author: Julien Dutheil
# Created: 12/11/2010

MACRO(TEST_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})

  SET(${lib_name}_NAMES ${lib_name} ${lib_name}.lib ${lib_name}.dll)
  FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})
  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} found here:")
    MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
    MESSAGE("   dynamic libraries: ${${lib_name}_LIBRARY}")
    MESSAGE(WARNING "Library ${lib_name} is already installed in the system tree. Test will be built against it. This may lead to unexpected results. You may want to do 'make install' before 'make test', or remove the installed version.")
  ELSE()
    SET(${lib_name}_LIBRARY "-L../src -lbpp-phyl")
    SET(${lib_name}_INCLUDE_DIR "../src/")
  ENDIF()
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(TEST_FIND_LIBRARY)

#Find the bpp-phyl library library:
TEST_FIND_LIBRARY(LIBS bpp-phyl Bpp/Phyl/Tree.h)

ADD_EXECUTABLE(test_tree test_tree.cpp)
TARGET_LINK_LIBRARIES(test_tree ${LIBS})
ADD_TEST(test_tree "test_tree")

ADD_EXECUTABLE(test_models test_models.cpp)
TARGET_LINK_LIBRARIES(test_models ${LIBS})
ADD_TEST(test_models "test_models")

ADD_EXECUTABLE(test_detailed_simulations test_detailed_simulations.cpp)
TARGET_LINK_LIBRARIES(test_detailed_simulations ${LIBS})
ADD_TEST(test_detailed_simulations "test_detailed_simulations")

ADD_EXECUTABLE(test_simulations test_simulations.cpp)
TARGET_LINK_LIBRARIES(test_simulations ${LIBS})
ADD_TEST(test_simulations "test_simulations")

ADD_EXECUTABLE(test_parsimony test_parsimony.cpp)
TARGET_LINK_LIBRARIES(test_parsimony ${LIBS})
ADD_TEST(test_parsimony "test_parsimony")

ADD_EXECUTABLE(test_likelihood test_likelihood.cpp)
TARGET_LINK_LIBRARIES(test_likelihood ${LIBS})
ADD_TEST(test_likelihood "test_likelihood")

ADD_EXECUTABLE(test_likelihood_nh test_likelihood_nh.cpp)
TARGET_LINK_LIBRARIES(test_likelihood_nh ${LIBS})
ADD_TEST(test_likelihood_nh "test_likelihood_nh")

ADD_EXECUTABLE(test_likelihood_clock test_likelihood_clock.cpp)
TARGET_LINK_LIBRARIES(test_likelihood_clock ${LIBS})
ADD_TEST(test_likelihood_clock "test_likelihood_clock")

ADD_EXECUTABLE(test_mapping test_mapping.cpp)
TARGET_LINK_LIBRARIES(test_mapping ${LIBS})
ADD_TEST(test_mapping "test_mapping")

ADD_EXECUTABLE(test_mapping_codon test_mapping_codon.cpp)
TARGET_LINK_LIBRARIES(test_mapping_codon ${LIBS})
ADD_TEST(test_mapping_codon "test_mapping_codon")

ADD_EXECUTABLE(test_nhx test_nhx.cpp)
TARGET_LINK_LIBRARIES(test_nhx ${LIBS})
ADD_TEST(test_nhx "test_nhx")

ADD_EXECUTABLE(test_bowker test_bowker.cpp)
TARGET_LINK_LIBRARIES(test_bowker ${LIBS})
ADD_TEST(test_bowker "test_bowker")

IF(UNIX)
  SET_PROPERTY(TEST test_detailed_simulations test_simulations test_parsimony test_models test_likelihood test_likelihood_nh test_likelihood_clock test_tree test_mapping test_mapping_codon test_nhx test_bowker PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:../src")
ENDIF()

IF(APPLE)
  SET_PROPERTY(TEST test_detailed_simulations test_simulations test_parsimony test_models test_likelihood test_likelihood_nh test_likelihood_clock test_tree test_mapping test_mapping_codon test_nhx test_bowker PROPERTY ENVIRONMENT "DYLD_LIBRARY_PATH=$ENV{DYLD_LIBRARY_PATH}:../src")
ENDIF()

IF(WIN32)
  SET(ENV{PATH} "$ENV{PATH};..\\src")
ENDIF()


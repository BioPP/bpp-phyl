# CMake script for bpp-phyl unit tests
# Authors:
#   Julien Dutheil
#   Francois Gindraud (2017)
# Created: 12/11/2010

# Add all tests.
# Any .cpp file in test/ is considered to be a test.
# It will be compiled as a standalone program (must contain a main()).
# A test is considered to succeed if it returns EXIT_SUCCESS (usually 0).

# Tests are linked to the the bpp-seq-shared target, which should point to the newly compiled library.
# CMake will then set the correct paths duing compilation (-L and -I).

file (GLOB test_cpp_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
foreach (test_cpp_file ${test_cpp_files})
  # Add each test (named as the filename without extension)
  get_filename_component (test_name ${test_cpp_file} NAME_WE)
  add_executable (${test_name} ${test_cpp_file})
  target_link_libraries (${test_name} bpp-phyl-shared ${BPP_LIBS_SHARED})
  add_test (
    NAME ${test_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${test_name}
    )
endforeach (test_cpp_file)

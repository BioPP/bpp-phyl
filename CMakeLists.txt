# SPDX-FileCopyrightText: The Bio++ Development Group
#
# SPDX-License-Identifier: CECILL-2.1

cmake_minimum_required(VERSION 3.10)
project(bpp-phyl CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compile options
add_compile_options(-Wall -Weffc++ -Wshadow -Wconversion -Wno-unknown-warning-option)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        RelWithDebInfo
        CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )
endif(NOT CMAKE_BUILD_TYPE)

if(NOT NO_DEP_CHECK)
    set(NO_DEP_CHECK
        FALSE
        CACHE BOOL
        "Disable dependencies check for building distribution only."
        FORCE
    )
endif(NOT NO_DEP_CHECK)

if(NO_DEP_CHECK)
    message(
        STATUS
        "Dependencies checking disabled. Only distribution can be built."
    )
else(NO_DEP_CHECK)
    #### Optional Debug

    option(USE_DEBUG "Enter Debug Mode" OFF)

    if(USE_DEBUG)
        message(STATUS "Enter Debug mode")
        add_compile_options(-DDEBUG)
    endif()

    # Libtool-like version number
    # CURRENT:REVISION:AGE => file.so.(C-A).A.R
    # current:  The most recent interface number that this library implements.
    # revision: The implementation number of the current interface.
    # age:      The difference between the newest and oldest interfaces that this
    #           library implements.
    # In other words, the library implements all the interface numbers in the
    # range from number current - age to current.
    set(${PROJECT_NAME}_VERSION_CURRENT "13")
    set(${PROJECT_NAME}_VERSION_REVISION "0")
    set(${PROJECT_NAME}_VERSION_AGE "0")

    #static linkage?
    if(NOT BUILD_STATIC)
        set(BUILD_STATIC
            FALSE
            CACHE BOOL
            "Compile and install static libraries."
            FORCE
        )
    endif()
    if(BUILD_STATIC)
        message(STATUS "Static libraries requested.")
    endif()

    # Effective version number computation
    math(
        EXPR
        ${PROJECT_NAME}_VERSION_MAJOR
        "${${PROJECT_NAME}_VERSION_CURRENT} - ${${PROJECT_NAME}_VERSION_AGE}"
    )
    set(${PROJECT_NAME}_VERSION_MINOR ${${PROJECT_NAME}_VERSION_AGE})
    set(${PROJECT_NAME}_VERSION_PATCH ${${PROJECT_NAME}_VERSION_REVISION})
    set(${PROJECT_NAME}_VERSION
        "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}"
    )

    set(PROJECT_VERSION ${${PROJECT_NAME}_VERSION})

    # Find dependencies (add install directory to search)
    if(CMAKE_INSTALL_PREFIX)
        set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
    endif(CMAKE_INSTALL_PREFIX)

    include(GNUInstallDirs)
    find_package(bpp-seq 14.0.0 REQUIRED)

    # CMake package
    set(cmake-package-location ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
    include(CMakePackageConfigHelpers)
    configure_package_config_file(
        package.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/package.cmake
        INSTALL_DESTINATION ${cmake-package-location}
    )
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/package-version.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )
    install(
        FILES ${CMAKE_CURRENT_BINARY_DIR}/package.cmake
        DESTINATION ${cmake-package-location}
        RENAME ${PROJECT_NAME}-config.cmake
    )
    install(
        FILES ${CMAKE_CURRENT_BINARY_DIR}/package-version.cmake
        DESTINATION ${cmake-package-location}
        RENAME ${PROJECT_NAME}-config-version.cmake
    )

    # Use eigen
    find_package(Eigen3 REQUIRED NO_MODULE)

    # Define the libraries
    add_subdirectory(src)

    # Doxygen
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        add_custom_target(
            apidoc
            cp Doxyfile ${CMAKE_BINARY_DIR}/Doxyfile-build
            COMMAND
                echo "OUTPUT_DIRECTORY=${CMAKE_BINARY_DIR}" >>
                ${CMAKE_BINARY_DIR}/Doxyfile-build
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile-build
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )
        add_custom_target(
            apidoc-stable
            cp Doxyfile ${CMAKE_BINARY_DIR}/Doxyfile-stable
            COMMAND
                echo "OUTPUT_DIRECTORY=${CMAKE_BINARY_DIR}" >>
                ${CMAKE_BINARY_DIR}/Doxyfile-stable
            COMMAND
                echo "HTML_HEADER=header.html" >>
                ${CMAKE_BINARY_DIR}/Doxyfile-stable
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile-stable
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )
    endif(DOXYGEN_FOUND)
endif(NO_DEP_CHECK)

# Packager
set(CPACK_PACKAGE_NAME "libbpp-phyl")
set(CPACK_PACKAGE_VENDOR "Bio++ Development Team")
set(CPACK_PACKAGE_VERSION "3.0.0")
set(CPACK_PACKAGE_VERSION_MAJOR "3")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Bio++ Phylogenetics library")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSES/CECILL-2.1.txt")
set(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
set(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
set(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
set(CPACK_SOURCE_IGNORE_FILES
    "/build/"
    "/html/"
    "/BppPhyl.tag"
    "/\\\\.git/"
    "/\\\\.gitignore"
    ${CPACK_SOURCE_IGNORE_FILES}
)

set(CPACK_SOURCE_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}"
)
set(CPACK_DEBSOURCE_PACKAGE_FILE_NAME
    "lib${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig"
)
include(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

if(NOT NO_DEP_CHECK)
    if(UNIX)
        #This creates rpm packages:
        add_custom_target(
            rpm
            rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz
        )
        add_dependencies(rpm dist)
    endif()

    set(CTEST_UPDATE_TYPE git)
    set(UPDATE_COMMAND "git")
    set(UPDATE_OPTIONS "")

    # Tests
    enable_testing()
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(test)
    endif(BUILD_TESTING)
endif(NOT NO_DEP_CHECK)

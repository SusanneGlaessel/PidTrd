cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(PidTrd)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_cpp_STANDARD 11)
set(PROJECT_VERSION 1.0)

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
endif()

set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "Using C++${CMAKE_CXX_STANDARD}")

message("CMAKE_PROJECT_NAME = ${CMAKE_PROJECT_NAME}")
message("PROJECT_NAME = ${PROJECT_NAME}")

include(ExternalProject)
message("ExternalProject = ${ExternalProject}")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)
#
# Options
#
#set(PIDTRD_BUNDLED_AT ON CACHE BOOL "Get and build AnalysisTree")
#set(PIDTRD_BUNDLED_AT_VERSION "master" CACHE STRING "Bundled AnalysisTree version")

#message("PIDTRD_BUNDLED_AT  ${PIDTRD_BUNDLED_AT}")
#message("PIDTRD_BUNDLED_AT_VERSION ${PIDTRD_BUNDLED_AT_VERSION}")

# Telling CMake where to find the ROOT installation.
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/lib")

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Core MathCore RIO Hist Tree)

# Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
add_definitions(${ROOT_CXX_FLAGS})

find_package(AnalysisTree)

message(STATUS "PROJECT_INCLUDE_DIRECTORIES ${PROJECT_INCLUDE_DIRECTORIES}")
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
message("ROOT_INCLUDE_DIRS = ${ROOT_INCLUDE_DIRS}")
set(PidTrd_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# by default build optimized code with debug symbols
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RELEASE)
endif ()

# in DEBUG mode make verbose Makefile
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif ()

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -DDEBUG -D__DEBUG -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ftree-vectorize -ffast-math -DNODEBUG")
message(STATUS "Using CXX flags for ${CMAKE_BUILD_TYPE}: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")

add_subdirectory(src)
if(AnalysisTree_FOUND)
    Message("AnalysisTree found. Corresponding interface will be
    compiled")
Else()
    Message("AnalysisTree is not found.")
EndIf()
add_subdirectory(interface)

 # Packaging routine
# For complete explanation
# @see https://cmake.org/cmake/help/git-master/manual/cmake-packages.7.html#creating-packages
include(GenerateExportHeader)

generate_export_header(PidTrd)
set_property(TARGET PidTrd PROPERTY VERSION ${PROJECT_VERSION})
set_property(TARGET PidTrd PROPERTY SOVERSION ${PROJECT_VERSION})
set_property(TARGET PidTrd PROPERTY INTERFACE_PidTrd_MAJOR_VERSION 1)
set_property(TARGET PidTrd APPEND PROPERTY COMPATIBLE_INTERFACE_STRING INTERFACE_PidTrd_MAJOR_VERSION)

install(TARGETS PidTrd EXPORT PidTrdTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        )

install(
        FILES
        ${HEADERS}
        "${CMAKE_CURRENT_BINARY_DIR}/pidtrd_export.h"
        DESTINATION
        include/pidtrd
        COMPONENT
        Devel
)

install(
        FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}_rdict.pcm"
        "${CMAKE_CURRENT_BINARY_DIR}/${PCM_FILE_NAME}.rootmap"
        DESTINATION
        lib
        OPTIONAL
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/PidTrd/PidTrdConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

export(EXPORT PidTrdTargets
        FILE "${CMAKE_CURRENT_BINARY_DIR}/PidTrd/PidTrdTargets.cmake")

message(STATUS "PidTrdTargets: ${PidTrd}")

set(ConfigPackageLocation lib/cmake/PidTrd)

set(PidTrd_INCLUDE_DIR "include")
set(PidTrd_LIBRARY_DIR "lib")

configure_package_config_file(
        PidTrdConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/PidTrd/PidTrdConfig.cmake
        INSTALL_DESTINATION ${ConfigPackageLocation}
        PATH_VARS PidTrd_INCLUDE_DIR PidTrd_LIBRARY_DIR
        #  [NO_SET_AND_CHECK_MACRO]
        #  [NO_CHECK_REQUIRED_COMPONENTS_MACRO]
        #  [INSTALL_PREFIX <path>]
)

install(EXPORT PidTrdTargets
        FILE
        PidTrdTargets.cmake
        DESTINATION
        ${ConfigPackageLocation}
        )

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/PidTrd/PidTrdConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/PidTrd/PidTrdConfigVersion.cmake"
        DESTINATION
        ${ConfigPackageLocation}
        COMPONENT
        Devel
        )

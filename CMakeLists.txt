cmake_minimum_required(VERSION 3.19)

project(most_common_string VERSION 0.0.1 LANGUAGES C CXX)

set(GTEST_GIT_REPO_PATH "${CMAKE_SOURCE_DIR}/3rd-party/gtest" CACHE STRING "Path for GTest git repository")

set(CMAKE_VERBOSE_MAKEFILE ON)

option(ENABLE_TEST "Build test" OFF)
option(ENABLE_BENCH "Build benchmarks" ON)
option(TARGET_ARCH "Target architecture" X86)
option(BUILD_STATIC "Build static" ON)
option(USE_OPENMP "Build with -qopenmp" ON)
option(BUILD_TYPE "Build type" Release)

include(cmake/common_configuration.cmake)

#
# Build library
#
add_subdirectory("${CMAKE_SOURCE_DIR}/lib")

#
# Build tests
#
if(ENABLE_TEST)
    add_subdirectory("${GTEST_GIT_REPO_PATH}")
    add_subdirectory("${CMAKE_SOURCE_DIR}/test")
endif()

#
# Build benchmarks
#
if(ENABLE_BENCH)
    add_subdirectory("${CMAKE_SOURCE_DIR}/perf")
endif()
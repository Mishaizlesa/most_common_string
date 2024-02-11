cmake_minimum_required(VERSION 3.19)

project(most_common_string VERSION 0.0.1 LANGUAGES C CXX)


set(CMAKE_VERBOSE_MAKEFILE ON)

option(ENABLE_TEST "Build test" OFF)
option(ENABLE_BENCH "Build benchmarks" ON)
option(TARGET_ARCH "Target architecture" RISCV_GENERIC)
option(BUILD_STATIC "Build static" ON)
option(USE_OPENMP "Build with -fopenmp" ON)
option(BUILD_TYPE "Build type" Release)

add_subdirectory(src)

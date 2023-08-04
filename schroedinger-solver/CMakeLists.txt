cmake_minimum_required(VERSION 3.16)

project(schroedinger-solver VERSION 0.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_CPPCHECK "cppcheck")

string(APPEND CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -fconcepts")
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer")
string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer")

include(ExternalProject)
include(FetchContent)

set(EIGEN_BUILD_DIR   ${CMAKE_BINARY_DIR}/eigen)
set(EIGEN_INSTALL_DIR ${CMAKE_SOURCE_DIR}/include/eigen3)
set(EIGEN_INCLUDE_DIR ${EIGEN_INSTALL_DIR})

set(SOURCES
    src/main.cpp
    src/solver_kernels.cpp
    src/solver_init_conditions.cpp
    src/potentials.cpp
    src/timer.cpp
)

ExternalProject_Add(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz
    PREFIX ${EIGEN_BUILD_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND
      ${CMAKE_COMMAND} -E copy_directory ${EIGEN_BUILD_DIR}/src/eigen/Eigen ${EIGEN_INCLUDE_DIR}/Eigen 
      && ${CMAKE_COMMAND} -E copy_directory ${EIGEN_BUILD_DIR}/src/eigen/unsupported ${EIGEN_INCLUDE_DIR}/unsupported
    TEST_COMMAND ""
)
include_directories(SYSTEM ${EIGEN_INCLUDE_DIR})

FetchContent_Declare(
  doctest
  GIT_REPOSITORY https://github.com/doctest/doctest.git
  GIT_TAG 7b9885133108ae301ddd16e2651320f54cafeba7
)
FetchContent_GetProperties(doctest)
if(NOT doctest_POPULATED)
  FetchContent_Populate(doctest)
endif()

find_package(OpenMP REQUIRED)

add_executable(schroedinger-solver ${SOURCES})
add_dependencies(schroedinger-solver eigen)
target_link_libraries(schroedinger-solver PRIVATE OpenMP::OpenMP_CXX)



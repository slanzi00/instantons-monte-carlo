cmake_minimum_required(VERSION 3.16)

project(lattice-solutions VERSION 0.1)

include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_CPPCHECK "cppcheck")
set(SOURCES
    correlators.cpp
    lattice.cpp
    main.cpp
    metropolis.cpp
)

string(APPEND CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -fconcepts")
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer")
string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer")

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

add_executable(lattice-solutions ${SOURCES})
target_link_libraries(lattice-solutions PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(lattice-solutions PRIVATE Boost::boost)

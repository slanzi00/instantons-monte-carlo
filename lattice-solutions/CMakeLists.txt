cmake_minimum_required(VERSION 3.16)

project(lattice-solutions VERSION 0.1)

include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_CPPCHECK "cppcheck")
set(SOURCES
    src/lattice.cpp
    src/correlators.cpp
    src/metropolis.cpp
    src/main.cpp
)

string(APPEND CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -fconcepts")
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer")
string(APPEND CMAKE_EXE_LINKER_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer")

find_package(Boost REQUIRED)
find_package(fmt REQUIRED)

add_executable(lattice-solutions ${SOURCES})
target_link_libraries(lattice-solutions PRIVATE Boost::boost)
target_link_libraries(lattice-solutions PRIVATE fmt::fmt)

add_custom_command(
    TARGET lattice-solutions POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_SOURCE_DIR}/plots
    COMMENT "Creating empty folders data and plots"
)

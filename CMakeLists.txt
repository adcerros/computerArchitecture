cmake_minimum_required(VERSION 3.10)
project(proy-grav-seq)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CPP_COMPILER g++)
add_definitions("-Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -O3 -DNDEBUG")

add_executable(sim-aos sim-aos.cpp)
add_executable(sim-soa sim-soa.cpp)

cmake_minimum_required(VERSION 3.26)
project(two_triangulations)

set(CMAKE_CXX_STANDARD 17)

add_executable(two_triangulations
        include/CutTriangulation.hpp
        tests/main_tests.cpp)

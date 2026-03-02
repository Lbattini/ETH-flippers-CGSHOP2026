cmake_minimum_required(VERSION 3.26)
project(2sat_linux)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OpenMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


add_executable(2sat_linux main.cpp)

set_target_properties(2sat_linux PROPERTIES
        CXX_CLANG_TIDY ""
)

target_link_libraries(2sat_linux PRIVATE OpenMP::OpenMP_CXX)

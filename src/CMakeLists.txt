cmake_minimum_required(VERSION 3.22)
project(A_star_for_flipdistance)

set(CMAKE_CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_package(CGAL REQUIRED COMPONENTS Core)

#CGAL_Qt6 is needed for the drawing.
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(centralTriangulation
        inOut.cpp
        korderDelaunay.cpp
        main.cpp
        geometryUtils.cpp
        genRandomTriangulations.cpp
)


#target_compile_definitions(korderDelaunay PUBLIC LARGEINTS=0)


target_include_directories(centralTriangulation
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
#        ${GUROBI_INCLUDE_DIRS}
)

target_link_libraries(centralTriangulation
        CGAL::CGAL
        CGAL::CGAL_Basic_viewer
)

#if(CGAL_Qt6_FOUND)
#
#    #link it with the required CGAL libraries
#    target_link_libraries(centralTriangulation PUBLIC CGAL::CGAL_Basic_viewer)
#endif()

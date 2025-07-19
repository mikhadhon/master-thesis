cmake_minimum_required(VERSION 3.31)
project(FlowComplex)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)

FetchContent_Declare(
        polyscope
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
        GIT_TAG        master
)

FetchContent_MakeAvailable(polyscope)

add_executable(
        FlowComplex main.cpp
        Delaunay.cpp
        Delaunay.h
        utils.cpp
        utils.h
)

find_package(CGAL CONFIG REQUIRED)
target_link_libraries(FlowComplex CGAL::CGAL polyscope)

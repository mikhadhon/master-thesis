cmake_minimum_required(VERSION 3.31)
project(FlowComplex4)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)

FetchContent_Declare(
        polyscope
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
        GIT_TAG        master
)

FetchContent_MakeAvailable(polyscope)

add_executable(
        FlowComplex4 main.cpp
        Delaunay.cpp
        Delaunay.h
        utils.cpp
        utils.h
        FlowComplex.cpp
        FlowComplex.h
)

find_package(CGAL CONFIG REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(FlowComplex4 CGAL::CGAL polyscope Eigen3::Eigen)

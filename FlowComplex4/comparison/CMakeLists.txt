cmake_minimum_required(VERSION 3.31)
project(Comparison)

set(CMAKE_CXX_STANDARD 20)

add_executable(
        Comparison main.cpp
)

find_package(CGAL CONFIG REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
target_link_libraries(Comparison CGAL::CGAL Eigen3::Eigen)

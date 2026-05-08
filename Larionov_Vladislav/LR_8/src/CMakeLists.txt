cmake_minimum_required(VERSION 3.14)
project(MyProject)

find_package(Matplot++ REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main Matplot++::matplot)

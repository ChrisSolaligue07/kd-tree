cmake_minimum_required(VERSION 3.23)
project(KdTree)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(CMAKE_CURRENT_SOURCE_DIR/src)

add_executable(kdtree_test src/main.cpp)
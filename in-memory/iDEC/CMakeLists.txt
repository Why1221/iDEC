cmake_minimum_required(VERSION 3.0)

project(idec)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-O0 -g")

include_directories(../Common)
add_executable(idec src/main.cc src/iDEC.cpp)

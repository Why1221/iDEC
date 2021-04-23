cmake_minimum_required(VERSION 3.0)

project(idec-external)

file(GLOB SRCS "*.cpp")

add_definitions(-DDEBUG)
add_executable(rt-example rtree-example.cc ${SRCS})
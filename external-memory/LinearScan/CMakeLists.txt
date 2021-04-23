cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 14)

project(linear-scan)

include_directories(../Common)


find_package(HighFive 2.1 REQUIRED)
add_executable(linear-scan-clion linear_scan_main.cc linear_scan_compact.cpp)
target_link_libraries(linear-scan-clion HighFive)
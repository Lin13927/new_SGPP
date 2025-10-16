cmake_minimum_required(VERSION 3.31)
project(SGPP)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_executable(SGPP main.cpp
        statistic.h)

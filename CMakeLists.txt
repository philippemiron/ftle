cmake_minimum_required(VERSION 3.2)
project(ftle)

if (APPLE)
    set(CMAKE_CXX_COMPILER /usr/local/bin/g++-5)
endif (APPLE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -lblas -llapack -fopenmp -std=c++11")

include_directories(include)
include_directories(src)

set(SOURCE_FILES
        include/exprtk.hpp
        main/ftle.cpp
        src/flowmap.cpp
        src/flowmap.h
        src/ftle.h
        src/rk45.h
        src/lyapunov.cpp
        src/lyapunov.h
        src/matrix.cpp
        src/matrix.h
        src/parameter.cpp
        src/parameter.h)

add_executable(ftle ${SOURCE_FILES})
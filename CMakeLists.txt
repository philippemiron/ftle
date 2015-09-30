cmake_minimum_required(VERSION 3.2)
project(ftle)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -lblas -llapack -fopenmp -std=c++11")

include_directories(include)
include_directories(src)

set(SOURCE_FILES
    include/exprtk.hpp
    main/ftle.cpp
    src/enum.h
    src/flowmap.cpp
    src/flowmap.h
    src/ftle.h
    src/lyapunov.cpp
    src/lyapunov.h
    src/matrix.cpp
    src/matrix.h
    src/parameter.cpp
    src/parameter.h
    src/writedata.cpp
    src/writedata.h)

add_executable(ftle ${SOURCE_FILES})
target_link_libraries(ftle ${CMAKE_SOURCE_DIR}/lib/libtecio_linux.a)
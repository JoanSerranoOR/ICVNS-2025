cmake_minimum_required(VERSION 3.29)
project(Grasp2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_C_FLAGS_DEBUG "-O0 -g")

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

add_executable(Grasp2
        main.cpp
        problem.cpp
        problem.h
        solution.cpp
        solution.h
        heuristic.cpp
        heuristic.h
        globals.cpp
        globals.h
)


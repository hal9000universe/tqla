cmake_minimum_required(VERSION 3.29)
project(tqla)

set(CMAKE_CXX_STANDARD 20)

add_executable(tqla main.cpp
        helper.h
        helper.cpp
        lossy_matrix.h
        lossy_matrix.cpp
        matrix.h
        matrix.cpp
        performance.h
        performance.cpp
        sparse_entry.h
        sparse_entry.cpp)

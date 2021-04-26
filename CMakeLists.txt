# File: CMakeLists
# 
# Project name:
# Contact:
# Author: 

cmake_minimum_required(VERSION 2.8.7 FATAL_ERROR)
project(reachsolver CXX C)

# Settings
set(PROJECT_FULLNAME "ReachSolver")
set(PROJECT_DESCRIPTION "")
string(TIMESTAMP reachsolver_MINORVERSION %m)
string(TIMESTAMP reachsolver_MAJORVERSION %y)
set(hypro_VERSION "${reachsolver_MAJORVERSION}.${reachsolver_MINORVERSION}")
set(PROJECT_VERSION "${reachsolver_VERSION}")
set(reachsolver_NAME "ReachSolver" )
set(reachsolver_DESCRIPTION ${PROJECT_DESCRIPTION})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_SOURCE_DIR}/cmake/Modules/")
include(macros)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Options
option(STATICLIB_SWITCH "Compile a statically linked version of the library." OFF)


# Libraries
load_library(${PROJECT_NAME} GTest 0)

enable_testing()

# Subdirectories
add_subdirectory(src)
add_subdirectory(src/test)

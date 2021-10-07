# File: CMakeLists
# 
# Project name:
# Contact:
# Author: 

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(reachsolver CXX C)
set(CMAKE_VERBOSE_MAKEFILE on)

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

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Options
option(ENABLE_EXAMPLES "Enable generating examples" ON)
option(ENABLE_TESTING "Enable Unit Testing" ON)
option(ENABLE_DOC "Enable generating documents" OFF)

# Libraries
load_library(${PROJECT_NAME} GTest 0)
find_package(autodiff)

enable_testing()

# Subdirectories
include_directories(src)
add_subdirectory(src)

if(ENABLE_EXAMPLES)
  message("Enable generating examples")
  add_subdirectory(examples)
endif()

if(ENABLE_TESTING)
  message("Enable unit testing")
  enable_testing()
  add_subdirectory(src/test)
endif()

if(ENABLE_DOC)
  message("Enable generating documents")
  add_subdirectory(doc)
endif()

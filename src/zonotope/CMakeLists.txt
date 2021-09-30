# Create object
file(GLOB_RECURSE sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
list(FILTER sources EXCLUDE REGEX old.cpp)



add_custom_target(allTests
        DEPENDS ${PROJECT_NAME}
)

set(test_sources
        zonotope_test.cpp
)

SET_SOURCE_FILES_PROPERTIES(${test_sources} PROPERTIES COMPILE_FLAGS -O0)

# Unit Tests

add_executable(runUnitTests ${test_sources})

target_include_directories(runUnitTests PRIVATE ${GTEST_INCLUDE_DIR})

target_link_libraries(
        runUnitTests LINK_PUBLIC
        ${GTEST_LIBRARIES}
        ${PROJECT_NAME}-shared
)

#cotire(runUnitTests)

add_test(NAME test COMMAND runUnitTests)

add_dependencies(allTests runUnitTests)
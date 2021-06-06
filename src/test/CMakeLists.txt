set(test_sources
        ZonotopeTest.cpp
)

SET_SOURCE_FILES_PROPERTIES(${test_sources} PROPERTIES COMPILE_FLAGS -O0)

# Unit Tests

add_executable(runUnitTests ${test_sources})

target_include_directories(runUnitTests PRIVATE ${GTEST_INCLUDE_DIR})

target_link_libraries(
        runUnitTests LINK_PUBLIC
        ${PROJECT_NAME}-shared
        gtest_main
        gtest
        pthread
)

add_test(NAME test COMMAND runUnitTests)

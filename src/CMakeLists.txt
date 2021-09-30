# Create static and shared libraries

file(GLOB_RECURSE sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
list(FILTER sources EXCLUDE REGEX old.cpp|ReachSpecification.c|test/*)

#add_subdirectory(zonotope)
#add_subdirectory(dynamics)

add_library(${PROJECT_NAME}-static STATIC ${sources})
set_target_properties( ${PROJECT_NAME}-static PROPERTIES
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION}"
    OUTPUT_NAME ${PROJECT_NAME}
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    LINKER_LANGUAGE CXX
)
target_link_libraries(${PROJECT_NAME}-static PUBLIC ${${PROJECT_NAME}_LIBRARIES_STATIC})
target_link_libraries(${PROJECT_NAME}-static PUBLIC pthread)

add_library(${PROJECT_NAME}-shared SHARED ${sources})
set_target_properties( ${PROJECT_NAME}-shared PROPERTIES
    VERSION "${PROJECT_VERSION}"
    SOVERSION "${PROJECT_VERSION}"
    OUTPUT_NAME ${PROJECT_NAME}
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    LINKER_LANGUAGE CXX
)
target_link_libraries(${PROJECT_NAME}-shared PUBLIC ${${PROJECT_NAME}_LIBRARIES_DYNAMIC})
target_link_libraries(${PROJECT_NAME}-shared PUBLIC pthread)

########################################################################################
#cotire(${PROJECT_NAME}-shared)
#cotire(${PROJECT_NAME}-static)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.cpp")
install(TARGETS ${PROJECT_NAME}-static
        EXPORT ${PROJECT_NAME}
        LIBRARY DESTINATION ${LIB_INSTALL_DIR}
        ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/static)
install(TARGETS ${PROJECT_NAME}-shared
        EXPORT ${PROJECT_NAME}
        LIBRARY DESTINATION ${LIB_INSTALL_DIR}
        ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/static)
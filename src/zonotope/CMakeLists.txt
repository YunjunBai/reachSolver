# Create object library first
file(GLOB_RECURSE sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
list(FILTER sources EXCLUDE REGEX old.cpp)

#add_library(${PROJECT_NAME}-objects OBJECT ${sources})
#target_include_dirs_from(${PROJECT_NAME}-objects SYSTEM PUBLIC EIGEN3 Boost_SHARED)
#target_include_directories(${PROJECT_NAME}-objects
#		PRIVATE ${${PROJECT_NAME}_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}
#		PUBLIC
#		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#		$<INSTALL_INTERFACE:include/${PROJECT_NAME}>
#	)

# Create static and shared libraries
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


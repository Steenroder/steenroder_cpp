cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(steenroder)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "./cmake")

# Project settings
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
set(EXECUTABLE_PATH "${PROJECT_SOURCE_DIR}/bin")

# Header file directories
include_directories("include")

# Google Test system Libraries
set(GTEST_ROOT ${CMAKE_SOURCE_DIR}/tests/lib/googletest)
if(NOT EXISTS ${GOOGLETEST_DIR})
    find_package(Git REQUIRED)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} clone https://github.com/google/googletest ${GTEST_ROOT}
    )
endif()

find_package(GTest)
if (NOT GTest_FOUND)
    message(WARNING
            "GTest was not found. Tests will not be available.")
    return()
endif ()

find_package(GTest REQUIRED)
include_directories("${GTEST_INCLUDE_DIR}")
#find_package(OpenMP REQUIRED)

# Doxygen
option(BUILD_DOC "Build Documentation" ON)

find_package(Doxygen)
if(DOXYGEN_FOUND)
    set(BUILD_DOC_DIR ${CMAKE_SOURCE_DIR}/build/docs)
    if(NOT EXISTS ${BUILD_DOC_DIR})
        file(MAKE_DIRECTORY ${BUILD_DOC_DIR})
    endif()

    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    message("Doxygen buld started")
    add_custom_target(Doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
else(DOXYGEN_FOUND)
    message("Doxygen needs to be installed to generate the documentation.")
endif(DOXYGEN_FOUND)

# Setup subdirectories
enable_testing()
add_subdirectory(docs)
add_subdirectory(tests)

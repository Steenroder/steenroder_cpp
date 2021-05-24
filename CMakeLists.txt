cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(steenroder)
set(CMAKE_CXX_STANDARD 11)

# Project settings
set(SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")
set(EXECUTABLE_PATH "${PROJECT_SOURCE_DIR}/bin")

# Header file directories
include_directories("include")
include_directories("external/AnyOption")

# Unit tests using Google Test
option(BUILD_TESTS "Build Unit Tests" ON)
if(BUILD_TESTS)
  enable_testing()
  include(GoogleTest)
  add_subdirectory(tests)
endif()

# Documentation using Doxygen
option(BUILD_DOCS "Build Documentation" ON)
if(BUILD_DOCS)
  add_subdirectory(docs)
endif()

# Third-party libraries
#find_package(OpenMP REQUIRED)

# Build Targets
add_subdirectory(src)

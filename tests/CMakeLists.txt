include(AddGoogleTest)
add_subdirectory("${PROJECT_SOURCE_DIR}/external/googletest" "external/googletest")

# Keeps cache clean
mark_as_advanced(
  BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
  gmock_build_tests gtest_build_samples gtest_build_tests
  gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
  )

# Keeps IDEs clean
set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)
set_target_properties(gmock PROPERTIES FOLDER extern)
set_target_properties(gmock_main PROPERTIES FOLDER extern)

# Test creation function
macro(steenroder_add_test test_name target_file)
  add_executable(${test_name} ${target_file})
  target_link_libraries(${test_name} gtest gmock gtest_main)
  gtest_discover_tests(${test_name}
    WORKING_DIRECTORY ${EXECUTABLE_PATH}
    PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
    )
  set_target_properties(${test_name} PROPERTIES FOLDER tests)
endmacro()

# Targets
steenroder_add_test(example TestExample.cpp)

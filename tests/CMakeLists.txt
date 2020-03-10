# Test creation function
function(steenroder_add_test target_name target_file)
    add_executable(${target_name} ${target_file})
    target_link_libraries(${target_name}
      "${GTEST_BOTH_LIBRARIES}")
    target_compile_features(looper PRIVATE cxx_range_for)
    add_test(NAME ${target_name}
      WORKING_DIRECTORY ${EXECUTABLE_PATH}
      COMMAND ${target_name})
endfunction()

# Targets
steenroder_add_test(example TestExample.cpp)

add_custom_target(stn COMMENT "Builds targets.")
function(stn_target DATATYPE COEFF)
  set(target_name stn_${DATATYPE}_${COEFF})
  add_executable(${target_name} ../external/AnyOption/anyoption.cpp barcodes.cpp)
  target_compile_definitions(${target_name} PRIVATE
    DATATYPE=${DATATYPE}
    COEFF="${COEFF}")
  add_dependencies(stn ${target_name})
endfunction()

stn_target(double 2)



function(stn_dualize DATATYPE COEFF)
  set(target_name dualize_${DATATYPE}_${COEFF})
  add_executable(${target_name} ../external/AnyOption/anyoption.cpp dualize.cpp)
  target_compile_definitions(${target_name} PRIVATE
    DATATYPE=${DATATYPE}
    COEFF="${COEFF}")
  add_dependencies(stn ${target_name})
endfunction()

stn_dualize(double 2)

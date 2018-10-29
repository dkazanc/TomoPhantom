project(tomophantomMatlab)


find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM MX_LIBRARY ENG_LIBRARY )


set(MATLAB_DEST_DIR "" CACHE PATH "Directory of the Matlab wrappers")
if (MATLAB_DEST_DIR)
set(MATLAB_DEST "${MATLAB_DEST_DIR}")
else() 
set(MATLAB_DEST "${CMAKE_INSTALL_PREFIX}/matlab")
endif()
message(STATUS "Matlab wrappers will be installed in " ${MATLAB_DEST})

file(GLOB CPU_MEX_FILES
    "${CMAKE_SOURCE_DIR}/Wrappers/MATLAB/src/*.c"   
)

#message("CPU_MEX_FILES " ${CPU_MEX_FILES})

list(LENGTH CPU_MEX_FILES num)


MATH(EXPR num "${num}-1")
#set(num "-1")
message("found ${num} files")

foreach(tgt RANGE 0 ${num})
  message("number " ${tgt})
  list(LENGTH CPU_MEX_FILES num2)
  message("the list is ${num2}")
  #list(GET CPU_TARGETS ${tgt} current_target)
  list(GET CPU_MEX_FILES ${tgt} current_file_name)
  get_filename_component(current_file ${current_file_name} NAME)
  string(REGEX MATCH "(.+).c" match ${current_file})
  if (NOT ${match} EQUAL "" )
  set (current_target ${CMAKE_MATCH_1})
  endif()
  message("matlab_add_mex target " ${current_file} " and " ${current_target})
  matlab_add_mex(
    NAME ${current_target}
    SRC 
      ${current_file_name} 
            
    LINK_TO ctomophantom ${Matlab_LIBRARIES}
    )
    
target_include_directories(${current_target}
   PUBLIC ${CMAKE_SOURCE_DIR}/Core
   ${MATLAB_INCLUDE_DIR})
   set_property(TARGET ${current_target} PROPERTY C_STANDARD 99)
   list(APPEND CPU_MEX_TARGETS ${current_target})
   INSTALL(TARGETS ${current_target} DESTINATION "${MATLAB_DEST}")
endforeach()
   
add_custom_target(MatlabWrapper DEPENDS ${CPU_MEX_TARGETS})



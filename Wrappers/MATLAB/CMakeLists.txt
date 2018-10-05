project(tomophantomMatlab)


find_package(Matlab REQUIRED COMPONENTS MAIN_PROGRAM MX_LIBRARY ENG_LIBRARY )


set(MATLAB_DEST_DIR "" CACHE PATH "Directory of the Matlab wrappers")
if (MATLAB_DEST_DIR)
set(MATLAB_DEST "${MATLAB_DEST_DIR}")
else() 
set(MATLAB_DEST "${CMAKE_INSTALL_PREFIX}/matlab")
endif()
message(STATUS "Matlab wrappers will be installed in " ${MATLAB_DEST})
#C:\Users\ofn77899\Documents\Projects\CCPi\GitHub\CCPi-FISTA_Reconstruction\Core\regularisers_CPU
# matlab_add_mex(
    # NAME CPU_ROF
    # SRC 
      # ${CMAKE_SOURCE_DIR}/Wrappers/Matlab/mex_compile/regularisers_CPU/ROF_TV.c
    # LINK_TO cilreg ${Matlab_LIBRARIES}
    # )
    
# target_include_directories(CPU_ROF 
   # PUBLIC ${CMAKE_SOURCE_DIR}/Core/regularisers_CPU
   # ${CMAKE_SOURCE_DIR}/Core/regularisers_GPU
   # ${CMAKE_SOURCE_DIR}/Core/inpainters_CPU
   # ${CMAKE_SOURCE_DIR}/Core/
   # ${MATLAB_INCLUDE_DIR})
   
   # matlab_add_mex(
    # NAME CPU_TNV
    # SRC 
      # ${CMAKE_SOURCE_DIR}/Wrappers/Matlab/mex_compile/regularisers_CPU/TNV.c 
    # LINK_TO cilreg ${Matlab_LIBRARIES}
    # )
    
# target_include_directories(CPU_TNV 
   # PUBLIC ${CMAKE_SOURCE_DIR}/Core/regularisers_CPU
   # ${CMAKE_SOURCE_DIR}/Core/regularisers_GPU
   # ${CMAKE_SOURCE_DIR}/Core/inpainters_CPU
   # ${CMAKE_SOURCE_DIR}/Core/
   # ${MATLAB_INCLUDE_DIR})
   
#set (CPU_MEX_FILES "regularisers_CPU/TNV.c;regularisers_CPU/ROF_TV.c")
#set (MEX_TARGETS "CPU_TNV;CPU_ROF")
#list(APPEND MEX_TARGETS "CPU_TNV")
#list(APPEND MEX_TARGETS "CPU_ROF")

file(GLOB CPU_MEX_FILES
    "${CMAKE_SOURCE_DIR}/Wrappers/MATLAB/src/*.c"
    #"${CMAKE_SOURCE_DIR}/Wrappers/Matlab/mex_compile/regularisers_GPU/*.c"
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
   PUBLIC ${CMAKE_SOURCE_DIR}/Core/regularisers_CPU
   ${CMAKE_SOURCE_DIR}/Core/regularisers_GPU
   ${CMAKE_SOURCE_DIR}/Core/inpainters_CPU
   ${CMAKE_SOURCE_DIR}/Core/
   ${MATLAB_INCLUDE_DIR})
   set_property(TARGET ${current_target} PROPERTY C_STANDARD 99)
   list(APPEND CPU_MEX_TARGETS ${current_target})
   INSTALL(TARGETS ${current_target} DESTINATION "${MATLAB_DEST}")
endforeach()
   
add_custom_target(MatlabWrapper DEPENDS ${CPU_MEX_TARGETS})



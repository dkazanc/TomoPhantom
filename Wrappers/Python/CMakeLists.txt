#   Copyright 2017 Edoardo Pasca
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
find_package(PythonInterp REQUIRED)
if (PYTHONINTERP_FOUND)
  message ("Current Python " ${PYTHON_VERSION_STRING} " found " ${PYTHON_EXECUTABLE})
endif()
if (CONDA_BUILD)
else()
set(PYTHON_DEST_DIR "" CACHE PATH "Directory of the Matlab wrappers")
  if (PYTHON_DEST_DIR)
   set(PYTHON_DEST "${PYTHON_DEST_DIR}")
  else() 
    set(PYTHON_DEST "${CMAKE_INSTALL_PREFIX}/python")
  endif()
  message(STATUS "Python wrappers will be installed in " ${PYTHON_DEST})
endif()

find_package(PythonInterp)
find_package(PythonLibs)
if (PYTHONINTERP_FOUND)
  message(STATUS "Found PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}")
  message(STATUS "Python version ${PYTHON_VERSION_STRING}")
endif()
if (PYTHONLIBS_FOUND)
  message(STATUS "Found PYTHON_INCLUDE_DIRS=${PYTHON_INCLUDE_DIRS}")
  message(STATUS "Found PYTHON_LIBRARIES=${PYTHON_LIBRARIES}")
endif()

if (PYTHONINTERP_FOUND)
    message("Python found " ${PYTHON_EXECUTABLE})
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    #set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/module/__init__.py")
    set (DEPS       "${CMAKE_BINARY_DIR}/Core/")
    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")
    
    if (WIN32)
      set (LIBRARY_LIB ${CMAKE_BINARY_DIR}/Core/${CMAKE_BUILD_TYPE})
      if (CONDA_BUILD)
	    set (LIBRARY_LIB ${CMAKE_BINARY_DIR}/Core/)
      endif()
    else()
      set (LIBRARY_LIB ${CMAKE_BINARY_DIR}/Core/)
    endif()

    message(STATUS "LIBRARY_LIB = ${LIBRARY_LIB}")
    configure_file(${SETUP_PY_IN} ${SETUP_PY})

        if (CONDA_BUILD)
      add_custom_command(OUTPUT ${OUTPUT}
                       COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src
                       COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/tomophantom ${CMAKE_CURRENT_BINARY_DIR}/tomophantom
                       COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/PhantomLibrary/models ${CMAKE_CURRENT_BINARY_DIR}/tomophantom
                       COMMAND ${CMAKE_COMMAND} -E env CIL_VERSION=${CIL_VERSION}
                                                       PREFIX=${CMAKE_SOURCE_DIR}/Core 
                                                       LIBRARY_INC=${CMAKE_SOURCE_DIR}/Core 
                                                       LIBRARY_LIB=${CMAKE_BINARY_DIR}/Core
                                                       ${PYTHON_EXECUTABLE} ${SETUP_PY} install
                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                       DEPENDS ctomophantom)
      add_custom_target(PythonWrapper ALL DEPENDS ${OUTPUT})

    else()
      add_custom_command(OUTPUT ${OUTPUT}
                       COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_BINARY_DIR}/src
                       COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/tomophantom ${CMAKE_CURRENT_BINARY_DIR}/tomophantom
                       COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/PhantomLibrary/models ${CMAKE_CURRENT_BINARY_DIR}/tomophantom
                       COMMAND ${CMAKE_COMMAND} -E env CIL_VERSION=${CIL_VERSION}
                                                       PREFIX=${CMAKE_SOURCE_DIR}/Core 
                                                       LIBRARY_INC=${CMAKE_SOURCE_DIR}/Core 
                                                       LIBRARY_LIB=${CMAKE_BINARY_DIR}/Core/${CMAKE_BUILD_TYPE}
                                                       ${PYTHON_EXECUTABLE} ${SETUP_PY} build_ext --inplace
                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                       DEPENDS ctomophantom
                       )
                       
      add_custom_target(PythonWrapper ALL DEPENDS ${OUTPUT})
      install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tomophantom 
              DESTINATION ${PYTHON_DEST})
      if (WIN32)
      install(FILES ${CMAKE_BINARY_DIR}/Core/${CMAKE_BUILD_TYPE}/ctomophantom.dll 
              DESTINATION ${PYTHON_DEST})
      else()
          install(FILES ${CMAKE_BINARY_DIR}/Core/libctomophantom.so
            DESTINATION ${PYTHON_DEST}
            CONFIGURATIONS ${CMAKE_BUILD_TYPE})
      endif()
      #install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/wip/test.py
      #              ${CMAKE_CURRENT_SOURCE_DIR}/wip/dvc_input.json
	#	    ${CMAKE_SOURCE_DIR}/Core/wip/central_grid.roi
    #          DESTINATION ${PYTHON_DEST})
    #  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/tomophantom/__init__.py
    #          DESTINATION ${PYTHON_DEST}/ccpi)
    endif()
    
    
    

    #install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install)")
endif()

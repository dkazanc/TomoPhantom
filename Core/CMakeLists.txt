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

project(tomophantom)

message("CIL_VERSION ${CIL_VERSION}")

find_package(Python 3.8 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

# Detect the installed nanobind package and import it into CMake
# Remove execute_process when scipy-build-core is used
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
   set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
   set (CMAKE_STATIC_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_STATIC_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
   
endif()

message("Building TomoPhantom core modules as a shared library")
message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")
message("CMAKE_C_FLAGS ${CMAKE_C_FLAGS}")
message("CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS}")
message("CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS}")
message("CMAKE_STATIC_LINKER_FLAGS ${CMAKE_STATIC_LINKER_FLAGS}")

message("Building TomoPhantom as a shared library")
nanobind_add_module(tomophantom "SHARED"
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DModel_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DModelSino_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP3DModel_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP3DModelSino_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DSinoNum_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/utils.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/bindings.c
		)
target_link_libraries(tomophantom ${EXTRA_LIBRARIES} )
include_directories(tomophantom PUBLIC 
                      ${LIBRARY_INC}/include 
           	      ${CMAKE_CURRENT_SOURCE_DIR}
		              )

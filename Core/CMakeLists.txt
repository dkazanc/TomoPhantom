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

set(CMAKE_BUILD_TYPE "Release")

if(WIN32)
  set (FLAGS "/DWIN32 /EHsc /DCCPi_EXPORTS /openmp:experimental")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
  set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:MSVCRT.lib")  
  set (EXTRA_LIBRARIES)		
  message("library lib: ${LIBRARY_LIB}")  
elseif(UNIX)
   set (FLAGS "-O2 -funsigned-char -Wall -Wl,--no-undefined  -DCCPiReconstructionIterative_EXPORTS")  
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
   set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
   if (OPENMP_FOUND)
     set (EXTRA_LIBRARIES 
		"gomp"
		"m"
		)
    else()
       set (EXTRA_LIBRARIES 
		"m"
		)
    endif()
endif()
message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")

message("Building TomoPhantom as a shared library")
add_library(tomophantom "SHARED"
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DModel_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DModelSino_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP3DModel_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP3DModelSino_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DSinoNum_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/utils.c
		)
target_link_libraries(tomophantom ${EXTRA_LIBRARIES} )
include_directories(tomophantom PUBLIC 
                      ${LIBRARY_INC}/include 
           	      ${CMAKE_CURRENT_SOURCE_DIR}
		              )

## Install
set(PYTHON_SITE_DIR ${CMAKE_INSTALL_PREFIX}/source)
set(CMAKE_POSITION_INDEPENDENT_CODE
    ON
    CACHE BOOL "Position independent code" FORCE)
#set(INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin)
#set(INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)
#set(INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)

if (UNIX)
message ("I'd install into ${CMAKE_INSTALL_PREFIX}/lib")
install(TARGETS tomophantom
	LIBRARY DESTINATION lib
	CONFIGURATIONS ${CMAKE_BUILD_TYPE} 
	)
elseif(WIN32)
message ("I'd install into ${CMAKE_INSTALL_PREFIX} lib bin")
  install(TARGETS tomophantom 
	RUNTIME DESTINATION bin
	ARCHIVE DESTINATION lib
	CONFIGURATIONS ${CMAKE_BUILD_TYPE} 
	)
endif()

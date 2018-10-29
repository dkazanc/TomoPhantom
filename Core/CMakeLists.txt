#   Copyright 2018 Edoardo Pasca
#cmake_minimum_required (VERSION 3.0)

project(RGL_core)
#https://stackoverflow.com/questions/13298504/using-cmake-with-setup-py

# The version number.

#set (CIL_VERSION $ENV{CIL_VERSION} CACHE INTERNAL "Core Imaging Library version" FORCE)

# conda orchestrated build
message("CIL_VERSION ${CIL_VERSION}")
#include (GenerateExportHeader)


find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
   set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
   set (CMAKE_STATIC_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_STATIC_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
   
endif()

## Build the regularisers package as a library
message("Building TomoPhantom as a shared library")

message("CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}")
message("CMAKE_C_FLAGS ${CMAKE_C_FLAGS}")
message("CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS}")
message("CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS}")
message("CMAKE_STATIC_LINKER_FLAGS ${CMAKE_STATIC_LINKER_FLAGS}")

set(CMAKE_BUILD_TYPE "Release")

if(WIN32)
  set (FLAGS "/DWIN32 /EHsc /DCCPi_EXPORTS")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
  set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:MSVCRT.lib")
  
  set (EXTRA_LIBRARIES)
		
  message("library lib: ${LIBRARY_LIB}")
  
elseif(UNIX)
   set (FLAGS "-O2 -funsigned-char -Wall  -Wl,--no-undefined  -DCCPiReconstructionIterative_EXPORTS ")  
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

## Build the regularisers package as a library
message("Building TomoPhantom as a shared library")

#set(CMAKE_C_COMPILER /apps/pgi/linux86-64/17.4/bin/pgcc)
#set(CMAKE_C_FLAGS "-acc -Minfo -ta=tesla:cc20 -openmp")
#set(CMAKE_C_FLAGS "-acc -Minfo -ta=multicore -openmp -fPIC")
add_library(ctomophantom "SHARED"
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DModel_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DModelSino_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP3DModel_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP3DModelSino_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DSinoNum_core.c
	    ${CMAKE_CURRENT_SOURCE_DIR}/utils.c
		)
#add_library(ctomophantom_static "STATIC"
#	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DModel_core.c
#	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DModelSino_core.c
#	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP3DModel_core.c
#	    ${CMAKE_CURRENT_SOURCE_DIR}/utils.c
#	    ${CMAKE_CURRENT_SOURCE_DIR}/TomoP2DSinoNum_core.c
#	    )
target_link_libraries(ctomophantom ${EXTRA_LIBRARIES} )
include_directories(ctomophantom PUBLIC 
                      ${LIBRARY_INC}/include 
					  ${CMAKE_CURRENT_SOURCE_DIR}
		              )

## Install

set(INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin)
set(INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)
set(INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)

if (UNIX)
message ("I'd install into ${CMAKE_INSTALL_PREFIX}/lib")
install(TARGETS ctomophantom
	LIBRARY DESTINATION lib
	CONFIGURATIONS ${CMAKE_BUILD_TYPE} 
	)
elseif(WIN32)
message ("I'd install into ${CMAKE_INSTALL_PREFIX} lib bin")
  install(TARGETS ctomophantom 
	RUNTIME DESTINATION bin
	ARCHIVE DESTINATION lib
	CONFIGURATIONS ${CMAKE_BUILD_TYPE} 
	)
endif()


#install(TARGETS ctomophantom
#	LIBRARY  DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
#	PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" COMPONENT dev
#	RUNTIME  DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
#	ARCHIVE  DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
#   CONFIGURATIONS ${CMAKE_BUILD_TYPE}
#	)

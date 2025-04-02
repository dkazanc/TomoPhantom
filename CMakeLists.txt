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

cmake_minimum_required(VERSION 3.15...3.27)

if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()

project(tomophantom)

# The version number.
set (CIL_VERSION_MAJOR 0)
set (CIL_VERSION_MINOR 11)
set (CIL_VERSION_PATCH 0)

set (CIL_VERSION '${CIL_VERSION_MAJOR}.${CIL_VERSION_MINOR}.${CIL_VERSION_PATCH}' CACHE INTERNAL "Core Imaging Library version" FORCE)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option (BUILD_MATLAB_WRAPPER "Build Matlab Wrappers" OFF)
option (BUILD_PYTHON_WRAPPER "Build Python Wrappers" OFF)
option (CONDA_BUILD "Conda Build" OFF)

add_subdirectory(Core)


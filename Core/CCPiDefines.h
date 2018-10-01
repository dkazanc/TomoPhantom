/**
# -*- coding: utf-8 -*-
#   This work is part of the Core Imaging Library developed by
#   Visual Analytics and Imaging System Group of the Science Technology
#   Facilities Council, STFC

#   Copyright 2018 CCPi

#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at

#       http://www.apache.org/licenses/LICENSE-2.0

#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

#   Code is derived from code developed by Prof. Brian Bay
*/

#ifndef CCPIDEFINES_H
#define CCPIDEFINES_H

#if defined(_WIN32) || defined(__WIN32__)
  #if defined(cilDVC_EXPORTS) ||  defined(CCPiNexusWidget_EXPORTS) || defined(CCPi_EXPORTS)//  add by CMake 
    #define  CCPI_EXPORT __declspec(dllexport)
    #define EXPIMP_TEMPLATE
  #else
    #define  CCPI_EXPORT __declspec(dllimport)
    #define EXPIMP_TEMPLATE extern
  #endif  CCPi_EXPORTS 
#elif defined(linux) || defined(__linux) || defined(__APPLE__)
 #define CCPI_EXPORT
#endif

//	Revised:

#define DAY_REV 25
#define MONTH_REV "Nov"
#define YEAR_REV 2017
	
#define VERSION 1.30

#endif

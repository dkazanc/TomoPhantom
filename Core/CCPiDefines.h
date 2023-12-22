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

#endif

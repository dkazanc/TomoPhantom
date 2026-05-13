#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

import ctypes
import ctypes.util
import sys, platform
import warnings


def c_shared_lib(lib_name, error=True):
    """Imports the library as a pre-built shared object"""
    if platform.system() == "Windows":
        load_dll = ctypes.windll.LoadLibrary
    else:
        load_dll = ctypes.cdll.LoadLibrary

    # Find the library and load it
    mylib_path = ctypes.util.find_library(lib_name)

    if mylib_path is not None:
        try:
            # No error if mylib_path is None; error if library name wrong
            return load_dll(mylib_path)
        except Exception as exc:
            if error:
                raise ModuleNotFoundError(lib_name) from exc
    warnings.warn(lib_name, ImportWarning, stacklevel=2)


from tomophantom.TomoP2D import *

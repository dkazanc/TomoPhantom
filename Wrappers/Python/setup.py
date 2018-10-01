#!/usr/bin/env python
"""
Copyright 2017  Srikanth Nagella / Daniil Kazantsev/ Edoardo Pasca

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

#import setuptools
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

#import os
import numpy
import platform	
#import sys

version = '1.0'
extra_include_dirs = [numpy.get_include(), '../functions/']
extra_library_dirs = []
extra_compile_args = []
extra_link_args    = []
extra_libraries    = []

if platform.system() == 'Windows':
    extra_compile_args += ['/DWIN32', '/openmp']
else:
    extra_compile_args += ['-fopenmp', '-O2', '-Wall', '-std=c99']
    extra_libraries += ['m','gomp']
    
setup(
    name='tomophantom',
    description='This is to generate phantom datasets for tomography experiments',
    version = version,
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize([ Extension("tomophantom.TomoP3D",
                            sources = [ "src/TomoP3D.pyx",
                                        "../functions/TomoP3DModel_core.c",
                                        "../functions/utils.c"
                                      ],
                            include_dirs = extra_include_dirs,
                            library_dirs = extra_library_dirs,
                            extra_compile_args = extra_compile_args,
                            libraries = extra_libraries,
                            extra_link_args = extra_link_args)]),
    zip_safe = False,
    include_package_data=True,
    package_data={'tomophantom':['*.dat','../functions/models/*.dat']},
    packages = {'tomophantom'}
)

setup(
    packages = {'tomophantom', 'tomophantom.supp'}
)


setup(
    name='tomophantom',
    description='This is to generate phantom datasets for tomography experiments',
    version = version,
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize([ Extension("tomophantom.TomoP2D",
                            sources = [ "src/TomoP2D.pyx",
                                        "../functions/TomoP2DModel_core.c",
                                        "../functions/TomoP2DModelSino_core.c",
                                        "../functions/utils.c"
                                      ],
                            include_dirs = extra_include_dirs,
                            library_dirs = extra_library_dirs,
                            extra_compile_args = extra_compile_args,
                            libraries = extra_libraries,
                            extra_link_args = extra_link_args)]),
    zip_safe = False,
    include_package_data=True,
    package_data={'tomophantom':['*.dat','../functions/models/*.dat']},
    packages = {'tomophantom'}
)

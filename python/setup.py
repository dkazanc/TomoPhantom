#!/usr/bin/env python

import setuptools
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import os
import numpy
import platform	
import sys

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
    extra_libraries += ['m']
    
setup(
    name='tomophantom',
    description='This is to generate phantom datasets for tomography experiments',
    version = version,
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize([ Extension("tomophantom.phantom3d",
                            sources = [ "src/phantom3d.pyx",
                                        "../functions/buildPhantom3D_core.c",
                                        "../functions/buildSino3D_core.c"
                                      ],
                            include_dirs = extra_include_dirs,
                            library_dirs = extra_library_dirs,
                            extra_compile_args = extra_compile_args,
                            libraries = extra_libraries,
                            extra_link_args = extra_link_args)]),
    zip_safe = False,
    include_package_data=True,
    package_data={'tomophantom':['*.dat','models/*.dat']},
    packages = {'tomophantom'}
)
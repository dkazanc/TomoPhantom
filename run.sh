#!/bin/bash  
echo "Building Tomophantom using CMake"  
# rm -r build
# Requires Cython, install it first: 
# pip install cython
mkdir build
cd build/
make clean
# install Python modules only
cmake ../ -DBUILD_PYTHON_WRAPPER=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install
# install Matlab modules only
# cmake ../ -DCONDA_BUILD=OFF -DBUILD_MATLAB_WRAPPER=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install
# install for both Matlab and Python modules only
# cmake ../ -DCONDA_BUILD=OFF -DBUILD_MATLAB_WRAPPER=ON -DBUILD_PYTHON_WRAPPER=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install
# in some cases if Matlab not found you need to give a path to ROOT dir, e.g.:
# cmake ../ -DCONDA_BUILD=OFF -DMatlab_ROOT_DIR=/home/algol/matlab2016/ -DBUILD_MATLAB_WRAPPER=ON -DBUILD_PYTHON_WRAPPER=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install
make install
cd install/python/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:.
# spyder
# one can run Matlab in Linux as:
# PATH="/path/to/mex/:$PATH" LD_LIBRARY_PATH="/path/to/library:$LD_LIBRARY_PATH" matlab

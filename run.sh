#!/bin/bash  
echo "Building Tomophantom using CMake"  
#rm -r build
mkdir build
cd build
#make clean
cmake ../ -DBUILD_PYTHON_WRAPPER=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=./install
make install
cd install/python/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:.
#spyder

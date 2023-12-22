.. _ref_installation:

Installation Guide
------------------
TomoPhantom consists of two parts: C-library that is build 
as a shared object using Cmake and Ctypes bindings and 
pure Python part that can be installed separately. 
Bellow we provide various recipes how software can be 
installed. We recommend to use TomoPhantom in Python as 
:ref:`ref_matlab` is not maintained at the moment.

.. _ref_python:

Python
======

Install TomoPhantom as a pre-built conda Python package:
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. code-block:: console

   $ conda install -c httomo tomophantom

or install with dependencies to perform reconstruction:

.. code-block:: console

   $ conda install -c httomo -c astra-toolbox tomophantom tomobar astra-toolbox ccpi-regularisation-cupy

Developers environment
+++++++++++++++++++++++
This sets the development environment to work in-place on the code.
It installs the library from the `httomo` conda channel.
If one needs to build the library first, please see how to 
:ref:`ref_libtomophantom`.

.. code-block:: console
    
   $ conda install -c httomo libtomophantom # get the library only
   $ git clone git@github.com/dkazanc/TomoPhantom.git # clone the repo   
   $ pip install -e .[dev] # the editable environment
   $ pytests tests/ # all tests should pass

Conda build
+++++++++++++
First one needs to build the `libtomophantom` C-library, the recipes for which are located in `conda-recipe_library/` folder:

.. code-block:: console
    
   $ export CIL_VERSION=3.0 # OR set CIL_VERSION=3.0 for Windows
   $ git clone git@github.com/dkazanc/TomoPhantom.git # clone the repo
   $ conda build conda-recipe_library # conda-build the library
   $ conda install path/to/the/tarball

Second, build the Python part of the library:

.. code-block:: console
    
   $ conda build conda-recipe -c httomo
   $ conda install path/to/the/tarball

.. _ref_libtomophantom:

Install libtomophantom from sources
+++++++++++++++++++++++++++++++++++++
With `Cmake <https://cmake.org>`_ one can build the `libtomophantom` library, which is a set of C functions. 
On Unix systems you will need a C compiler (gcc) and VC compiler on Windows. Note that in our 
`httomo <https://anaconda.org/httomo/>`_ conda channel we provide the pre-built binaries of libtomophantom.

Here an example of the `libtomophantom` library build on Linux. The shared library object 
will be installed into the conda lib location, 
see `echo $CONDA_PREFIX`.

.. code-block:: console
    
   $ git clone https://github.com/dkazanc/TomoPhantom.git
   $ mkdir build
   $ cd build
   $ cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
   $ cmake --build . # on Windows add a flag --config Release
   $ cmake --install .

The library should now be installed in your conda environment, you can `pip install` the Python part. 

.. _ref_matlab:

Matlab
======
Initially TomoPhantom was created as a Matlab application with mex compiled Matlab wrappers.
After 2018, however, the developers decided to stop supporting Matlab. We believe it still can be 
compiled and installed in Matlab, but its functionality is not garanteed. We recommend using it in :ref:`ref_python`. 
For Matlab installations please see the older `versions <https://github.com/dkazanc/TomoPhantom/releases>`_ of TomoPhantom.



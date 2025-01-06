.. _ref_installation:

Installation Guide
------------------

TomoPhantom consists of two parts: C-library that is build 
as a shared object using Cmake and Ctypes bindings and 
pure Python part that can be installed separately.

Bellow, we provide various recipes on how software can be 
installed. We recommend to use TomoPhantom in Python as 
:ref:`ref_matlab` is not maintained at the moment.

.. _ref_python:

Python
======

Install TomoPhantom as a pre-built conda Python package:
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This should be the first choice of quickly installing TomoPhantom in your *new* conda environment: 

.. code-block:: console

   $ conda install httomo::tomophantom

One can also install TomoPhantom with the reconstruction tools as dependencies:

.. code-block:: console

   $ $ conda install -c httomo -c conda-forge tomophantom tomobar astra-toolbox ccpi-regulariser

The above should work for Linux and Windows, however, the Mac OS users need to build TomoPhantom following this :ref:`ref_installation_condabuild` guide.


Developers environment
+++++++++++++++++++++++

This sets the development environment to work in-place on the **Python part** of code.
It installs the :code:`libtomophantom` library from the `httomo` conda channel. If one needs to build the library first, please see how to 
:ref:`ref_libtomophantom`.

.. code-block:: console
    
   $ conda install httomo::libtomophantom # get the library only
   $ git clone git@github.com/dkazanc/TomoPhantom.git # clone the repo   
   $ pip install -e .[dev] # the editable environment
   $ pytests tests/ # all tests should pass

.. _ref_installation_condabuild:

Conda build
+++++++++++

First, one needs to build the :code:`libtomophantom` C-library, the recipes for which are located in the dedicated :code:`conda-recipe_library` folder:

.. code-block:: console
    
   $ export CIL_VERSION=3.1 # OR set CIL_VERSION=3.1 for Windows
   $ git clone git@github.com/dkazanc/TomoPhantom.git # clone the repo
   $ conda build conda-recipe_library # conda-build the libtomophantom library
   $ conda install path/to/the/tarball

Second, build the Python part of the library:

.. code-block:: console
    
   $ conda build conda-recipe
   $ conda install path/to/the/tarball
   $ pytests tests/ # all tests should pass

.. _ref_libtomophantom:

Install libtomophantom from sources
++++++++++++++++++++++++++++++++++++

With `Cmake <https://cmake.org>`_ one can build the `libtomophantom` library, which is a set of C-modules. 
On Unix systems you will need a C compiler (gcc) and VC compiler on Windows.

Here an example of the `libtomophantom` library build on Linux. The shared library object 
will be installed into the conda lib location, check it with `echo $CONDA_PREFIX`.

.. code-block:: console
    
   $ git clone https://github.com/dkazanc/TomoPhantom.git
   $ mkdir build
   $ cd build
   $ cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
   $ cmake --build . # on Windows add a flag --config Release
   $ cmake --install .

The library should now be installed in your conda environment, you can just `pip install` the Python part. 

.. _ref_matlab:

Matlab
======

 In the beginning, TomoPhantom was created mainly for Matlab users with mex-compiled Matlab wrappers.
 After 2018, however, the developers decided to stop supporting Matlab. We believe it still can be 
 compiled and installed in Matlab, but its functionality is not guaranteed. We recommend using TomoPhantom in :ref:`ref_python`. 
 For Matlab installations, please see the older `versions <https://github.com/dkazanc/TomoPhantom/releases>`_ of TomoPhantom.



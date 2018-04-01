"""
Cython recipe to create an interface to C-functions (3D version)

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
# cython and ctypes
import cython
import ctypes

# import numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

from enum import Enum

# declare the interface to the C code
cdef extern float TomoP3DModel_core(float *A, int ModelSelected, int N, char* ModelParametersFilename)
cdef extern float TomoP3DObject_core(float *A, int N, char *Object, float C0, float x0, float y0, float z0, float a, float b, float c, float psi1, float psi2, float psi3, int tt)
cdef extern float checkParams3D(int *params_switch, int ModelSelected, char *ModelParametersFilename)
#cdef extern float buildSino3D_core(float *A, int ModelSelected, int N, int P, float *Th, int AngTot, int CenTypeIn, char* ModelParametersFilename)
#cdef extern float buildSino3D_core_single(float *A, int N, int P, float *Th, int AngTot, int CenTypeIn, int Object, float C0, float x0, float y0, float z0, float a, float b, float c, float phi_rot)
    
cdef packed struct object_3d:
    char[22] Obj
    np.float32_t C0
    np.float32_t x0
    np.float32_t y0
    np.float32_t z0
    np.float32_t a
    np.float32_t b
    np.float32_t c
    np.float32_t psi1
    np.float32_t psi2
    np.float32_t psi3

class Objects3D(Enum):
    '''Enumeration with the available objects for 3D phantoms'''
    GAUSSIAN  = 'gaussian'
    PARABOLOID  = 'paraboloid'
    ELLIPSOID = 'ellipsoid'
    CONE   = 'cone'
    CUBOID      = 'cuboid'
    ELLIPCYLINDER = 'ellipticalcylinder'

@cython.boundscheck(False)
@cython.wraparound(False)
def Model(int model_id, int phantom_size, str model_parameters_filename):
    """
    Create stationary 3D model : Model(model_id, phantom_size,model_parameters_filename)
    
    Takes in a input model_id and phantom_size and returns a phantom-model (3D) of phantom_size x phantom_size x phantom_size of type float32 numpy array.
    
    param: model_parameters_filename -- filename for the model parameters
    param: model_id -- a model id from the functions file
    param: phantom_size -- a phantom size in each dimension.
    
    returns: numpy float32 phantom array    
    """
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] phantom = np.zeros([phantom_size, phantom_size, phantom_size], dtype='float32')
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[int, ndim=1, mode="c"] params
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    checkParams3D(&params[0], model_id, c_string)
    testParams3D(params) # check parameters and terminate before running the core
    if params[3] == 1:
        ret_val = TomoP3DModel_core(&phantom[0,0,0], model_id, phantom_size, c_string)
    else:
        print("The selected model is temporal (4D), use 'ModelTemporal' function instead")
    return phantom
@cython.boundscheck(False)
@cython.wraparound(False)
def ModelTemporal(int model_id, int phantom_size, str model_parameters_filename):
    """
    Create temporal 4D (3D + time) model :Model(model_id, phantom_size,model_parameters_filename)
    
    Takes in a input model_id and phantom_size and returns a phantom-model (4D) of phantom_size x phantom_size x phantom_size x Time-frames of type float32 numpy array.
    
    param: model_parameters_filename -- filename for the model parameters
    param: model_id -- a model id from the functions file
    param: phantom_size -- a phantom size in each dimension.
    
    returns: numpy float32 phantom array    
    """
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[int, ndim=1, mode="c"] params
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    checkParams3D(&params[0], model_id, c_string)
    testParams3D(params) # check parameters and terminate before running the core
    cdef np.ndarray[np.float32_t, ndim=4, mode="c"] phantom = np.zeros([params[3], phantom_size, phantom_size, phantom_size], dtype='float32')
    if params[3] == 1:
        print("The selected model is stationary (3D), use 'Model' function instead")
    else:
        ret_val = TomoP3DModel_core(&phantom[0,0,0,0], model_id, phantom_size, c_string)
    return phantom
@cython.boundscheck(False)
@cython.wraparound(False)
def Object(int phantom_size, objlist):
    """
    Object (phantom_size,objlist)
    
    Takes in a input object description (list) and phantom_size and returns a phantom-object (3D) of phantom_size x phantom_size x phantom_size of type float32 numpy array.
    
    param: obj_params -- object parameters list
    param: phantom_size -- a phantom size in each dimension.
    
    returns: numpy float32 phantom array    
    """
#    cdef Py_ssize_t i
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] phantom = np.zeros([phantom_size, phantom_size, phantom_size], dtype='float32')
    cdef float ret_val
    if type(objlist) is dict:
        objlist = [objlist]
        
    for obj in objlist:
        if testParamsPY(obj):
            
            objectName = bytes(obj['Obj'].value, 'ascii')
            
            ret_val = TomoP3DObject_core(&phantom[0,0,0], phantom_size,
                                        objectName, 
                                        obj['C0'], 
                                        obj['x0'], 
                                        obj['y0'], 
                                        obj['z0'], 
                                        obj['a'], 
                                        obj['b'], 
                                        obj['c'], 
                                        obj['phi1'], 
                                        obj['phi2'], 
                                        obj['phi3'], 0)
    return phantom

def testParamsPY(obj):
    '''Performs a simple type check of the input parameters and a range check'''
    if not type(obj) is dict:
        raise TypeError('obj is not a dict {0}'.format(type(obj)))
    # type check
    typecheck = type(obj['C0']) is float
    if not typecheck:
        raise TypeError('C0 is not a float')
    typecheck = typecheck and type(obj['x0']) is float
    if not typecheck:
        raise TypeError('x0 is not a float')
    typecheck = typecheck and type(obj['y0']) is float
    if not typecheck:
        raise TypeError('y0 is not a float')
    typecheck = typecheck and type(obj['z0']) is float
    if not typecheck:
        raise TypeError('z0 is not a float')
    typecheck = typecheck and type(obj['a']) is float
    if not typecheck:
        raise TypeError('a is not a float')
    typecheck = typecheck and type(obj['b']) is float
    if not typecheck:
        raise TypeError('b is not a float')
    typecheck = typecheck and type(obj['c']) is float
    if not typecheck:
        raise TypeError('c is not a float')
    typecheck = typecheck and type(obj['phi1']) is float
    if not typecheck:
        raise TypeError('phi1 is not a float')
    typecheck = typecheck and type(obj['phi2']) is float
    if not typecheck:
        raise TypeError('phi2 is not a float')
    typecheck = typecheck and type(obj['phi3']) is float
    if not typecheck:
        raise TypeError('phi3 is not a float')
    
    # range check    
    rangecheck = obj['x0'] >= -1 and obj['x0'] <= 1
    if not rangecheck:
        raise ValueError('x0 is out of range. Must be between -1 and 1')
    rangecheck = rangecheck and obj['y0'] >= -1 and obj['y0'] <= 1
    if not rangecheck:
        raise ValueError('y0 is out of range. Must be between -1 and 1')
    rangecheck = rangecheck and obj['z0'] >= -1 and obj['z0'] <= 1
    if not rangecheck:
        raise ValueError('z0 is out of range. Must be between -1 and 1')
    rangecheck = rangecheck and obj['a'] > 0 and obj['a'] <= 2
    if not rangecheck:
        raise ValueError('a (object size) must be positive in [0,2] range')
    rangecheck = rangecheck and obj['b'] > 0 and obj['b'] <= 2
    if not rangecheck:
        raise ValueError('b (object size) must be positive in [0,2] range')
    rangecheck = rangecheck and obj['c'] > 0 and obj['c'] <= 2
    if not rangecheck:
        raise ValueError('c (object size) must be positive in [0,2] range')
    return rangecheck and typecheck

def testParams3D(obj):
    if obj[0] == 0:
        raise TypeError('Check if the library file <Phantom3DLibrary.dat> exists, the given path is correct and the syntax is valid')
    if obj[1] == 0:
        raise TypeError('The given model is not found, check available models in <Phantom3DLibrary.dat> file')
    if obj[2] == 0:
        raise TypeError('Components number cannot be negative, check <Phantom3DLibrary.dat> file')
    if obj[3] == 0:
        raise TypeError('TimeSteps cannot be negative, check <Phantom3DLibrary.dat> file')
    if obj[4] == 0:
        raise TypeError('Unknown name of the object, check <Phantom3DLibrary.dat> file')
    if obj[5] == 0:
        raise TypeError('C0 should not be equal to zero, check <Phantom3DLibrary.dat> file')
    if obj[6] == 0:
        raise TypeError('x0 (object position) must be in [-1,1] range, check <Phantom3DLibrary.dat> file')
    if obj[7] == 0:
        raise TypeError('y0 (object position) must be in [-1,1] range, check <Phantom3DLibrary.dat> file')
    if obj[8] == 0:
        raise TypeError('z0 (object position) must be in [-1,1] range, check <Phantom3DLibrary.dat> file')
    if obj[9] == 0:
        raise TypeError('a (object size) must be positive in [0,2] range, check <Phantom3DLibrary.dat> file')
    if obj[10] == 0:
        raise TypeError('b (object size) must be positive in [0,2] range, check <Phantom3DLibrary.dat> file')
    if obj[11] == 0:
        raise TypeError('c (object size) must be positive in [0,2] range, check <Phantom3DLibrary.dat> file')
    return 0
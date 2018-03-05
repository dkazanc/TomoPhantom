"""
Cython recipe to create an interface to C-functions (2D version)

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

# declare the interface to the C code
cdef extern float TomoP2DModel_core(float *A, int ModelSelected, int N, char* ModelParametersFilename, int platform)
cdef extern float TomoP2DObject(float *A, int N, char *Object, float C0, float x0, float y0, float a, float b, float phi_rot)
cdef extern float TomoP2DModelSino_core(float *A, int ModelSelected, int N, int P, float *Th, int AngTot, int CenTypeIn, char* ModelParametersFilename, int platform)
cdef extern float TomoP2DObjectSino(float *A, int N, int P, float *Th, int AngTot, int CenTypeIn, char *Object, float C0, float x0, float y0, float a, float b, float phi_rot)
cdef extern float extractSteps(int *steps, int ModelSelected, char *ModelParametersFilename)

cdef packed struct object_2d:
    char[16] Obj
    np.float32_t C0
    np.float32_t x0
    np.float32_t y0
    np.float32_t a
    np.float32_t b
    np.float32_t phi_rot
    
@cython.boundscheck(False)
@cython.wraparound(False)
def Model(int model_id, int phantom_size, str model_parameters_filename):
    """
    To generate stationary (2D) Model(model_id, phantom_size,model_parameters_filename)
    
    Takes in a input model_id and phantom_size and returns a phantom-model of phantom_size x phantom_size of type float32 numpy array.
    
    param: model_parameters_filename -- filename for the model parameters
    param: model_id -- a model id from the functions file
    param: phantom_size -- a phantom size in each dimension.
    
    returns: numpy float32 phantom array
    """
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] phantom = np.zeros([phantom_size, phantom_size], dtype='float32')
    cdef np.ndarray[int, ndim=1, mode="c"] steps
    steps = np.ascontiguousarray(np.zeros([1], dtype=ctypes.c_int))
    extractSteps(&steps[0], model_id, c_string)
    if steps[0] == 1:
        ret_val = TomoP2DModel_core(&phantom[0,0], model_id, phantom_size, c_string, 1)
    else:
        print("The selected model is temporal (3D), use 'ModelTemporal' function instead")
    return phantom
@cython.boundscheck(False)
@cython.wraparound(False)
def ModelTemporal(int model_id, int phantom_size, str model_parameters_filename):
    """
    to generate temporal (2D + time) Model(model_id, phantom_size,model_parameters_filename)
    
    Takes in a input model_id and phantom_size and returns a phantom-model of Time Frames x phantom_size x phantom_size of type float32 numpy array.
    
    param: model_parameters_filename -- filename for the model parameters
    param: model_id -- a model id from the functions file
    param: phantom_size -- a phantom size in each dimension.
    
    returns: numpy float32 phantom array
    """
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[int, ndim=1, mode="c"] steps
    steps = np.ascontiguousarray(np.zeros([1], dtype=ctypes.c_int))
    extractSteps(&steps[0], model_id, c_string)
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] phantom = np.zeros([steps[0], phantom_size, phantom_size], dtype='float32')    
    if steps[0] == 1:
        print("The selected model is static (2D), use 'Model' function instead")
    else:
        ret_val = TomoP2DModel_core(&phantom[0,0,0], model_id, phantom_size, c_string, 1)
    return phantom
@cython.boundscheck(False)
@cython.wraparound(False)
def Object(int phantom_size, object_2d[:] obj_params):
    """
    Object (phantom_size,object_parameters)
    
    Takes in a input object description (list) and phantom_size and returns a phantom-object of phantom_size x phantom_size of type float32 numpy array.
    
    param: phantom_size -- a phantom size in each dimension.
    param: obj_params -- object parameters list
    
    returns: numpy float32 phantom array
    
    """
    cdef Py_ssize_t i
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] phantom = np.zeros([phantom_size, phantom_size], dtype='float32')
    cdef float ret_val
    for i in range(obj_params.shape[0]):
        ret_val = TomoP2DObject(&phantom[0,0], phantom_size, obj_params[i].Obj, obj_params[i].C0, obj_params[i].x0, obj_params[i].y0, obj_params[i].a, obj_params[i].b, obj_params[i].phi_rot)
    return phantom
    
@cython.boundscheck(False)
@cython.wraparound(False)
def ModelSino(int model_id, int image_size, int detector_size, np.ndarray[np.float32_t, ndim=1, mode="c"] angles, str model_parameters_filename):
    
    """
    function to build a sinogram of the model (2D): ModelSino (model_id, image_size, detector_size, angles, model_parameters_filename)
    
    Takes in as input model_id, image_size, detector_size and projection angles and return a 2D sinogram corresponding to the model id.
    
    param: model_parameters_filename -- filename for the model parameters
    param: model_id -- a model id from the functions file
    param: image_size -- a phantom size in each dimension.
    param: detector_size -- int detector size.
    param: angles -- a numpy array of float values with angles in radians
    param: CenTypeIn -- 1 as default [0: radon, 1:astra]
    returns: numpy float32 phantom sinograms array.    
    np.flipud(np.fliplr(sinogram.transpose()))
    """
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] sinogram = np.zeros([detector_size,angles.shape[0]], dtype='float32')
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string    
    cdef int AngTot = angles.shape[0]
    cdef int CenTypeIn = 1 # astra center positioning
    cdef np.ndarray[int, ndim=1, mode="c"] steps
    steps = np.ascontiguousarray(np.zeros([1], dtype=ctypes.c_int))
    extractSteps(&steps[0], model_id, c_string)
    if steps[0] == 1:
        ret_val = TomoP2DModelSino_core(&sinogram[0,0], model_id, image_size, detector_size, &angles[0], AngTot, CenTypeIn, c_string, 1)
    else:
        print("The selected model is temporal (3D), use 'ModelSinoTemporal' function instead")
    return sinogram.transpose()
@cython.boundscheck(False)
@cython.wraparound(False)
def ModelSinoTemporal(int model_id, int image_size, int detector_size, np.ndarray[np.float32_t, ndim=1, mode="c"] angles, str model_parameters_filename):
    
    """
    function to build a 3D (2D +t) sinogram of the model: ModelSino (model_id, image_size, detector_size, angles, model_parameters_filename)
    
    Takes in as input model_id, image_size, detector_size and projection angles and return a 2D sinogram corresponding to the model id.
    
    param: model_parameters_filename -- filename for the model parameters
    param: model_id -- a model id from the functions file
    param: image_size -- a phantom size in each dimension.
    param: detector_size -- int detector size.
    param: angles -- a numpy array of float values with angles in radians
    param: CenTypeIn -- 1 as default [0: radon, 1:astra]
    returns: numpy float32 phantom sinograms array.    
    """    
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string    
    cdef int AngTot = angles.shape[0]
    cdef int CenTypeIn = 1 # astra center positioning
    cdef np.ndarray[int, ndim=1, mode="c"] steps
    steps = np.ascontiguousarray(np.zeros([1], dtype=ctypes.c_int))
    extractSteps(&steps[0], model_id, c_string)
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] sinogram = np.zeros([steps[0], detector_size, angles.shape[0]], dtype='float32')
    if steps[0] == 1:
        print("The selected model is stationary (2D), use 'ModelSino' function instead")
    else:
        ret_val = TomoP2DModelSino_core(&sinogram[0,0,0], model_id, image_size, detector_size, &angles[0], AngTot, CenTypeIn, c_string, 1)
    return sinogram
@cython.boundscheck(False)
@cython.wraparound(False)
def ObjectSino(int image_size, int detector_size, np.ndarray[np.float32_t, ndim=1, mode="c"] angles, object_2d[:] obj_params):    
    """
    ObjectSino (image_size, detector_size, angles,object parameters)
    
    Takes in as input object parameters list, image_size, detector_size and projection angles and return a 2D sinogram corresponding to the object
    
    param: image_size -- a phantom size in each dimension.
    param: detector_size -- int detector size.
    param: angles -- a numpy array of float values with angles in radians
    param: CenTypeIn -- 1 as default [0: radon, 1:astra]
    param: obj_params -- object parameters list
    returns: numpy float32 phantom sinograms array.
    
    """
    cdef Py_ssize_t i    
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] sinogram = np.zeros([detector_size,angles.shape[0]], dtype='float32')
    cdef float ret_val 
    cdef int AngTot = angles.shape[0]
    cdef int CenTypeIn = 1 # astra center posit
    for i in range(obj_params.shape[0]):
        ret_val = TomoP2DObjectSino(&sinogram[0,0], image_size, detector_size, &angles[0], AngTot, CenTypeIn, obj_params[i].Obj, obj_params[i].C0, obj_params[i].x0, obj_params[i].y0, obj_params[i].a, obj_params[i].b, obj_params[i].phi_rot)
    return sinogram.transpose()
"""
@cython.boundscheck(False)
@cython.wraparound(False)
def Object2(int phantom_size, paramlist):
    cdef Py_ssize_t i
    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] phantom = np.zeros([phantom_size, phantom_size], dtype='float32')
    cdef float ret_val
    for obj in obj_params:
        if testParams(obj):
            ret_val = TomoP2DObject(&phantom[0,0], phantom_size, 
                                    obj['Obj'], 
                                    obj['C0'], 
                                    obj['x0'], 
                                    obj['y0'], 
                                    obj['a'], 
                                    obj['b'], 
                                    obj['phi'])
    return phantom

def testParams(obj):
    typecheck = type(obj['Obj']) is str
    if not typecheck:
        raise TypeError('Obj is not a string')
    else:
        if not obj['Obj'] in ['gaussian', 'parabola', 'parabola1', 'ellipse', 'cone', 'rectangle']:
            raise ValueError('Model unknown: {0}'.format(obj['Obj']))
    typecheck = typecheck and type(obj['C0']) is float:
        if not typecheck:
            raise TypeError('C0 is not a float')
    typecheck = typecheck and type(obj['x0']) is float:
        if not typecheck:
            raise TypeError('x0 is not a float')
    typecheck = typecheck and type(obj['y0']) is float:
        if not typecheck:
            raise TypeError('y0 is not a float')
    typecheck = typecheck and type(obj['a']) is float:
        if not typecheck:
            raise TypeError('a is not a float')
    typecheck = typecheck and type(obj['b']) is float:
        if not typecheck:
            raise TypeError('b is not a float')
    typecheck = typecheck and type(obj['phi']) is float:
        if not typecheck:
            raise TypeError('phi is not a float')
           
    rangecheck = obj['x0'] >= -1 and obj['x0'] <= 1
    if not rangecheck:
        raise ValueError('x0 is out of range. Must be between -1 and 1')
    rangecheck = rangecheck and obj['y0'] >= -1 and obj['y0'] <= 1
    if not rangecheck:
        raise ValueError('y0 is out of range. Must be between -1 and 1')
    rangecheck = rangecheck and obj['a'] > 0
    if not rangecheck:
        raise ValueError('a is not positive.')
    rangecheck = rangecheck and obj['b'] > 0
    if not rangecheck:
        raise ValueError('b is not positive.')
    
    return rangecheck and typecheck
"""
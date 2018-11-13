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

from numbers import Number
from enum import Enum

import sys

# declare the interface to the C code
cdef extern float TomoP3DModel_core(float *A, int ModelSelected, long N1, long N2, long N3, long Z1, long Z2, char* ModelParametersFilename)
cdef extern float TomoP3DObject_core(float *A, long N1, long N2, long N3, long Z1, long Z2, char *Object, float C0, float x0, float y0, float z0, float a, float b, float c, float psi1, float psi2, float psi3, int tt)
cdef extern float TomoP3DModelSino_core(float *A, int ModelSelected, long Horiz_det, long Vert_det, long Z1, long Z2, long N, float *Theta_proj, int AngTot, char* ModelParametersFilename)
cdef extern float TomoP3DObjectSino_core(float *A, long Horiz_det, long Vert_det, long Z1, long Z2, long N, float *Angl_vector, int AngTot, char *Object,float C0, float x0, float y0, float z0, float a, float b, float c, float psi1, float psi2, float psi3, int tt)
cdef extern float checkParams3D(int *params_switch, int ModelSelected, char *ModelParametersFilename)
    
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
def Model(int model_id, phantom_size, str model_parameters_filename):
    """       
    Create stationary 3D model : Model(model_id, phantom_size, model_parameters_filename)
    
    Takes in an input model_id and phantom_size as a scalar or a tuple and returns a phantom-model (3D) of phantom_size of type float32 numpy array.
    
    param: model_parameters_filename -- filename for the model parameters
    param: model_id -- a model id from the functions file
    param: phantom_size -- a  scalar or a tuple with phantom dimesnsions. Can be phantom_size[1] (a scalar for the cubic phantom), or phantom_size[3] = [N1,N2,N3]  
       
    returns: numpy float32 phantom array    
    """
    cdef long N1,N2,N3
    if type(phantom_size) == tuple:
       N1,N2,N3 = [int(i) for i in phantom_size]
    else:
       N1 = N2 = N3 = phantom_size
    
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] phantom = np.zeros([N1, N2, N3], dtype='float32')
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[int, ndim=1, mode="c"] params
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    checkParams3D(&params[0], model_id, c_string)
    testParams3D(params) # check parameters and terminate before running the core
    if params[3] == 1:
        ret_val = TomoP3DModel_core(&phantom[0,0,0], model_id, N3, N2, N1, 0l, N1, c_string)
    else:
        print("The selected model is temporal (4D), use 'ModelTemporal' function instead")
    return phantom
@cython.boundscheck(False)
@cython.wraparound(False)
def ModelSub(int model_id, phantom_size, sub_index, str model_parameters_filename):
    """       
    Create a subset (vertical cutoff) of stationary 3D model : Model(model_id, phantom_size, sub_index, model_parameters_filename)
    
    Takes in an input model_id and phantom_size as a scalar or a tuple, then a sub_index a tuple defines selected vertical subset
    and returns a phantom-model (3D) of phantom_size of type float32 numpy array.    
    
    param: model_id -- a model id from the functions file
    param: sub_index -- a tuple containing 2 indeces [lower, upper] to select a vertical subset needed to be extracted
    param: phantom_size -- a  scalar or a tuple with phantom dimesnsions. Can be phantom_size[1] (a scalar for the cubic phantom), or phantom_size[3] = [N1,N2,N3]  
    param: model_parameters_filename -- filename for the model parameters
    
    returns: numpy float32 phantom array    
    """
    cdef long N1,N2,N3,Z1,Z2
    if type(phantom_size) == tuple:
       N1,N2,N3 = [int(i) for i in phantom_size]
    else:
       N1 = N2 = N3 = phantom_size
    
    Z1,Z2 = [int(i) for i in sub_index]
    
    rangecheck = Z2 > Z1
    if not rangecheck:
        raise ValueError('Upper index must be larger than the lower one')
    rangecheck = Z1 >= 0 and Z1 < N3
    if not rangecheck:
        raise ValueError('Range of the lower index is incorrect')
    rangecheck = Z2 >= 0 and Z2 <= N3
    if not rangecheck:
        raise ValueError('Range of the higher index is incorrect')
    
    sub_size = Z2 - Z1 # the size of the vertical slab
    
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] phantom = np.zeros([sub_size, N2, N3], dtype='float32')
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[int, ndim=1, mode="c"] params
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    checkParams3D(&params[0], model_id, c_string)
    testParams3D(params) # check parameters and terminate before running the core
    if params[3] == 1:
        ret_val = TomoP3DModel_core(&phantom[0,0,0], model_id, N3, N2, N1, Z1, Z2, c_string)
    else:
        print("The selected model is temporal (4D), use 'ModelTemporal' function instead")
    return phantom
@cython.boundscheck(False)
@cython.wraparound(False)
def ModelTemporal(int model_id, phantom_size, str model_parameters_filename):
    """
    Create temporal 4D (3D + time) model :Model(model_id, phantom_size,model_parameters_filename)
    
    Takes in a input model_id and phantom_size and returns a phantom-model (4D) of tuple size x Time-frames of type float32 numpy array.
    
    param: model_parameters_filename -- filename for the model parameters
    param: model_id -- a model id from the functions file
    param: phantom_size -- a  scalar or a tuple with phantom dimesnsions. Can be phantom_size[1] (a scalar for the cubic phantom), or phantom_size[3] = [N1,N2,N3]  
    
    returns: numpy float32 phantom array    
    """
    cdef long N1,N2,N3
    if type(phantom_size) == tuple:
       N1,N2,N3 = [int(i) for i in phantom_size]
    else:
       N1 = N2 = N3 = phantom_size
    
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[int, ndim=1, mode="c"] params
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    checkParams3D(&params[0], model_id, c_string)
    testParams3D(params) # check parameters and terminate before running the core
    cdef np.ndarray[np.float32_t, ndim=4, mode="c"] phantom = np.zeros([params[3], N1, N2, N3], dtype='float32')
    if params[3] == 1:
        print("The selected model is stationary (3D), use 'Model' function instead")
    else:
        ret_val = TomoP3DModel_core(&phantom[0,0,0,0], model_id, N3, N2, N1, 0l, N1, c_string)
    return phantom
@cython.boundscheck(False)
@cython.wraparound(False)
def ModelTemporalSub(int model_id, phantom_size, sub_index, str model_parameters_filename):
    """
    Create a subset of temporal 4D (3D + time) model: Model(model_id, phantom_size, sub_index, model_parameters_filename)
    
    Takes in a input model_id. phantom_size, sub_index a tuple defines selected vertical subset and returns a phantom-model (4D) of tuple size x Time-frames of type float32 numpy array.
    
    param: model_parameters_filename -- filename for the model parameters
    param: phantom_size -- a  scalar or a tuple with phantom dimesnsions. Can be phantom_size[1] (a scalar for the cubic phantom), or phantom_size[3] = [N1,N2,N3]  
    param: sub_index -- a tuple containing 2 indeces [lower, upper] to select a vertical subset needed to be extracted
    param: model_id -- a model id from the functions file 
    
    returns: numpy float32 phantom array
    """
    cdef long N1,N2,N3
    if type(phantom_size) == tuple:
       N1,N2,N3 = [int(i) for i in phantom_size]
    else:
       N1 = N2 = N3 = phantom_size
       
    Z1,Z2 = [int(i) for i in sub_index]
    
    rangecheck = Z2 > Z1
    if not rangecheck:
        raise ValueError('Upper index must be larger than the lower one')
    rangecheck = Z1 >= 0 and Z1 < N3
    if not rangecheck:
        raise ValueError('Range of the lower index is incorrect')
    rangecheck = Z2 >= 0 and Z2 <= N3
    if not rangecheck:
        raise ValueError('Range of the higher index is incorrect')
    
    sub_size = Z2 - Z1 # the size of the vertical slab
        
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[int, ndim=1, mode="c"] params
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    checkParams3D(&params[0], model_id, c_string)
    testParams3D(params) # check parameters and terminate before running the core
    cdef np.ndarray[np.float32_t, ndim=4, mode="c"] phantom = np.zeros([params[3], sub_size, N2, N3], dtype='float32')
    if params[3] == 1:
        print("The selected model is stationary (3D), use 'Model' function instead")
    else:
        ret_val = TomoP3DModel_core(&phantom[0,0,0,0], model_id, N3, N2, N1, Z1, Z2, c_string)
    return phantom
@cython.boundscheck(False)
@cython.wraparound(False)
def Object(phantom_size, objlist):
    """
    Object (phantom_size,objlist)
    
    Takes in a input object description (list) and phantom_size and returns a phantom-object (3D) of a tuple size of type float32 numpy array.
    
    param: obj_params -- object parameters list
    param: phantom_size -- a  scalar or a tuple with phantom dimesnsions. Can be phantom_size[1] (a scalar for the cubic phantom), or phantom_size[3] = [N1,N2,N3]  
    
    returns: numpy float32 phantom array    
    """
#    cdef Py_ssize_t i
    cdef long N1,N2,N3
    if type(phantom_size) == tuple:
       N1,N2,N3 = [int(i) for i in phantom_size]
    else:
       N1 = N2 = N3 = phantom_size

    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] phantom = np.zeros([N1, N2, N3], dtype='float32')
    cdef float ret_val
    if type(objlist) is dict:
        objlist = [objlist]
        
    for obj in objlist:
        if testParamsPY(obj):
            if sys.version_info.major == 3:
                objectName = bytes(obj['Obj'].value, 'ascii')
            elif sys.version_info.major == 2:
                objectName = bytes(obj['Obj'].value)

            ret_val = TomoP3DObject_core(&phantom[0,0,0], N3, N2, N1, 0l, N1,
                                        objectName, 
                                        float(obj['C0']),
                                        float(obj['y0']),
                                        float(obj['x0']),
                                        float(obj['z0']),
                                        float(obj['a']),
                                        float(obj['b']),
                                        float(obj['c']), 
                                        float(obj['phi1']),
                                        float(0.0), 
                                        float(0.0), 0)
    return phantom
@cython.boundscheck(False)
@cython.wraparound(False)
def ModelSino(int model_id, phantom_size, Horiz_det, Vert_det, np.ndarray[np.float32_t, ndim=1, mode="c"] angles, str model_parameters_filename):
    """  
    Creates 3D analytical projection data of the dimension [AngTot, Vert_det, Horiz_det] 
    
    Takes in a input model_id. phantom_size, detector sizes, array of projection angles
    
    param: model_id -- a model id from the functions file 
    param: phantom_size -- a  scalar or a tuple with phantom dimesnsions. Can be phantom_size[1] (a scalar for the cubic phantom)
    param: Horiz_det -- horizontal detector size
    param: Vert_det -- vertical detector size (currently Vert_det = phantom_size)
    params: angles -- 1D array of projection angles in degrees
    param: model_parameters_filename -- filename for the model parameters
    
    returns: numpy float32 phantom array (3D)
    """
    if type(phantom_size) == tuple:
       raise ValueError('Please give a scalar for phantom size, projection data cannot be obtained from non-cubic phantom')
    cdef int AngTot = angles.shape[0]
    
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] projdata = np.zeros([AngTot, Vert_det, Horiz_det], dtype='float32')
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[int, ndim=1, mode="c"] params
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    checkParams3D(&params[0], model_id, c_string)
    testParams3D(params) # check parameters and terminate before running the core
    if params[3] == 1:
        ret_val = TomoP3DModelSino_core(&projdata[0,0,0], model_id, Horiz_det, Vert_det, 0l, Vert_det, phantom_size, &angles[0], AngTot, c_string)
    else:
        print("The selected model is temporal (4D), use 'ModelTemporalSino' function instead")
    return np.swapaxes(projdata,0,1)
@cython.boundscheck(False)
@cython.wraparound(False)
def ModelSinoTemporal(int model_id, phantom_size, Horiz_det, Vert_det, np.ndarray[np.float32_t, ndim=1, mode="c"] angles, str model_parameters_filename):
    """  
    Creates 4D (3D + time) analytical projection data of the dimension [TimeFrames, AngTot, Vert_det, Horiz_det]
    
    Takes in a input model_id. phantom_size, detector sizes, array of projection angles
    
    param: model_id -- a model id from the functions file 
    param: phantom_size -- a  scalar or a tuple with phantom dimesnsions. Can be phantom_size[1] (a scalar for the cubic phantom)
    param: Horiz_det -- horizontal detector size
    param: Vert_det -- vertical detector size (currently Vert_det = phantom_size)
    params: angles -- 1D array of projection angles in degrees
    param: model_parameters_filename -- filename for the model parameters
    
    returns: numpy float32 phantom array (4D)
    """
    if type(phantom_size) == tuple:
       raise ValueError('Please give a scalar for phantom size, projection data cannot be obtained from non-cubic phantom')
    cdef int AngTot = angles.shape[0]
    
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[int, ndim=1, mode="c"] params
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    checkParams3D(&params[0], model_id, c_string)
    testParams3D(params) # check parameters and terminate before running the core
    cdef np.ndarray[np.float32_t, ndim=4, mode="c"] projdata = np.zeros([params[3], AngTot, Vert_det, Horiz_det], dtype='float32')
    if params[3] == 1:
        print("The selected model is stationary (3D), use 'ModelSino' function instead")
    else:
        ret_val = TomoP3DModelSino_core(&projdata[0,0,0,0], model_id, Horiz_det, Vert_det, 0l, Vert_det, phantom_size, &angles[0], AngTot, c_string)
    return np.swapaxes(projdata,1,2)
@cython.boundscheck(False)
@cython.wraparound(False)
def ModelSinoSub(int model_id, phantom_size, Horiz_det, Vert_det, sub_index, np.ndarray[np.float32_t, ndim=1, mode="c"] angles, str model_parameters_filename):
    """  
    Creates 3D analytical projection data of the dimension [AngTot, Subset_size, Horiz_det] 
    
    sub_index - defines a subset of the vertical detector dimension
    
    Takes in a input model_id. phantom_size, detector sizes, array of projection angles
    
    param: model_id -- a model id from the functions file 
    param: phantom_size -- a  scalar or a tuple with phantom dimesnsions. Can be phantom_size[1] (a scalar for the cubic phantom)
    param: Horiz_det -- horizontal detector size
    param: Vert_det -- vertical detector size (currently Vert_det = phantom_size)
    params: angles -- 1D array of projection angles in degrees
    param: model_parameters_filename -- filename for the model parameters
    
    returns: numpy float32 phantom array (3D)
    """
    if type(phantom_size) == tuple:
       raise ValueError('Please give a scalar for phantom size, projection data cannot be obtained from non-cubic phantom')
    cdef int AngTot = angles.shape[0]
    
    Z1,Z2 = [int(i) for i in sub_index]
    
    rangecheck = Z2 > Z1
    if not rangecheck:
        raise ValueError('Upper index must be larger than the lower one')
    rangecheck = Z1 >= 0 and Z1 < Vert_det
    if not rangecheck:
        raise ValueError('Range of the lower index is incorrect')
    rangecheck = Z2 >= 0 and Z2 <= Vert_det
    if not rangecheck:
        raise ValueError('Range of the higher index is incorrect')
    
    sub_size = Z2 - Z1 # the size of the vertical slab
    
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] projdata = np.zeros([AngTot, sub_size, Horiz_det], dtype='float32')
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[int, ndim=1, mode="c"] params
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    checkParams3D(&params[0], model_id, c_string)
    testParams3D(params) # check parameters and terminate before running the core
    if params[3] == 1:
        ret_val = TomoP3DModelSino_core(&projdata[0,0,0], model_id, Horiz_det, Vert_det, Z1, Z2, phantom_size, &angles[0], AngTot, c_string)
    else:
        print("The selected model is temporal (4D), use 'ModelTemporalSino' function instead")
    return np.swapaxes(projdata,0,1)
@cython.boundscheck(False)
@cython.wraparound(False)
def ModelSinoTemporalSub(int model_id, phantom_size, Horiz_det, Vert_det, sub_index, np.ndarray[np.float32_t, ndim=1, mode="c"] angles, str model_parameters_filename):
    """  
    Creates 4D analytical projection data of the dimension [AngTot, Subset_size, Horiz_det] 
    
    sub_index - defines a subset of the vertical detector dimension
    
    Takes in a input model_id. phantom_size, detector sizes, array of projection angles
    
    param: model_id -- a model id from the functions file 
    param: phantom_size -- a  scalar or a tuple with phantom dimesnsions. Can be phantom_size[1] (a scalar for the cubic phantom)
    param: Horiz_det -- horizontal detector size
    param: Vert_det -- vertical detector size (currently Vert_det = phantom_size)
    params: angles -- 1D array of projection angles in degrees
    param: model_parameters_filename -- filename for the model parameters
    
    returns: numpy float32 phantom array (3D)
    """
    if type(phantom_size) == tuple:
       raise ValueError('Please give a scalar for phantom size, projection data cannot be obtained from non-cubic phantom')
    cdef int AngTot = angles.shape[0]
    
    Z1,Z2 = [int(i) for i in sub_index]
    
    rangecheck = Z2 > Z1
    if not rangecheck:
        raise ValueError('Upper index must be larger than the lower one')
    rangecheck = Z1 >= 0 and Z1 < Vert_det
    if not rangecheck:
        raise ValueError('Range of the lower index is incorrect')
    rangecheck = Z2 >= 0 and Z2 <= Vert_det
    if not rangecheck:
        raise ValueError('Range of the higher index is incorrect')
    
    sub_size = Z2 - Z1 # the size of the vertical slab
    
    cdef float ret_val
    py_byte_string = model_parameters_filename.encode('UTF-8')
    cdef char* c_string = py_byte_string
    cdef np.ndarray[int, ndim=1, mode="c"] params
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    checkParams3D(&params[0], model_id, c_string)
    testParams3D(params) # check parameters and terminate before running the core
    cdef np.ndarray[np.float32_t, ndim=4, mode="c"] projdata = np.zeros([params[3], AngTot, sub_size, Horiz_det], dtype='float32')
    if params[3] == 1:
        print("The selected model is stationary (3D), use 'ModelSino' function instead")
    else:
        ret_val = TomoP3DModelSino_core(&projdata[0,0,0,0], model_id, Horiz_det, Vert_det, Z1, Z2, phantom_size, &angles[0], AngTot, c_string)
    return np.swapaxes(projdata,1,2)
@cython.boundscheck(False)
@cython.wraparound(False)
def ObjectSino(phantom_size, Horiz_det, Vert_det, np.ndarray[np.float32_t, ndim=1, mode="c"] angles, objlist):
    """  
    Creates 3D analytical projection data of the dimension [AngTot, Vert_det, Horiz_det] 
    
    Takes in a input model_id. phantom_size, detector sizes, array of projection angles

    param: phantom_size -- a  scalar or a tuple with phantom dimesnsions. Can be phantom_size[1] (a scalar for the cubic phantom)
    param: Horiz_det -- horizontal detector size
    param: Vert_det -- vertical detector size (currently Vert_det = phantom_size)
    params: angles -- 1D array of projection angles in degrees
    param: objlist -- a list of parameters for the object
    
    returns: numpy float32 phantom array (3D)
    """
    if type(phantom_size) == tuple:
       raise ValueError('Please give a scalar for phantom size, projection data cannot be obtained from non-cubic phantom')
    cdef int AngTot = angles.shape[0]
    
    cdef np.ndarray[np.float32_t, ndim=3, mode="c"] projdata = np.zeros([AngTot, Vert_det, Horiz_det], dtype='float32')
    cdef float ret_val
    if type(objlist) is dict:
        objlist = [objlist]
        
    for obj in objlist:
        if testParamsPY(obj):
            if sys.version_info.major == 3:
                objectName = bytes(obj['Obj'].value, 'ascii')
            elif sys.version_info.major == 2:
                objectName = bytes(obj['Obj'].value)
            
            if (("gaussian" == objectName) or ("paraboloid" == objectName) or ("ellipsoid" == objectName)):
                ret_val = TomoP3DObjectSino_core(&projdata[0,0,0], Horiz_det, Vert_det, 0l, Vert_det, phantom_size, &angles[0], AngTot,
                                        objectName, 
                                        float(obj['C0']),
                                        float(obj['y0']),
                                        float(-obj['z0']),
                                        float(-obj['x0']),
                                        float(obj['b']),
                                        float(obj['a']),
                                        float(obj['c']), 
                                        float(0.0),
                                        float(0.0), 
                                        float(obj['phi1']), 0)
            elif ("elliptical_cylinder" == objectName):
                ret_val = TomoP3DObjectSino_core(&projdata[0,0,0], Horiz_det, Vert_det, 0l, Vert_det, phantom_size, &angles[0], AngTot,
                                        objectName, 
                                        float(obj['C0']),
                                        float(obj['x0']),
                                        float(-obj['y0']),
                                        float(obj['z0']),
                                        float(obj['b']),
                                        float(obj['a']),
                                        float(obj['c']), 
                                        float(0.0),
                                        float(0.0), 
                                        float(obj['phi1']), 0)
            else:
                ret_val = TomoP3DObjectSino_core(&projdata[0,0,0], Horiz_det, Vert_det, 0l, Vert_det, phantom_size, &angles[0], AngTot,
                                        objectName, 
                                        float(obj['C0']),
                                        float(obj['x0']),
                                        float(obj['y0']),
                                        float(obj['z0']),
                                        float(obj['a']),
                                        float(obj['b']),
                                        float(obj['c']), 
                                        float(0.0),
                                        float(0.0), 
                                        float(-obj['phi1']), 0)
    return np.swapaxes(projdata,0,1)
@cython.boundscheck(False)
@cython.wraparound(False)

def testParamsPY(obj):
    '''Performs a simple type check of the input parameters and a range check'''
    if not type(obj) is dict:
        raise TypeError('obj is not a dict {0}'.format(type(obj)))
    # type check
    for k,v in obj.items():
        if not isinstance(v, Number):
            if not k == 'Obj':
                raise TypeError(k, 'is not a Number')
    typecheck = True
    
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

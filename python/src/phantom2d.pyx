"""
tomophantom.pyx
Phantom 2D Generator
"""

import cython

# import numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern float buildPhantom2D_core(float *A, int ModelSelected, int N, char* ModelParametersFilename)
cdef extern float buildPhantom2D_core_single(float *A, int N, int Object, float C0, float x0, float y0, float a, float b, float phi_rot)
#cdef extern float buildSino2D_core(float *A, int ModelSelected, int N, int P, float *Th, int AngTot, int CenTypeIn, char* ModelParametersFilename)
#cdef extern float buildSino2D_core_single(float *A, int N, int P, float *Th, int AngTot, int CenTypeIn, int Object, float C0, float x0, float y0, float a, float b, float phi_rot)
	
cdef packed struct object_2d:
	np.int_t Obj
	np.float32_t C0
	np.float32_t x0
	np.float32_t y0
	np.float32_t a
	np.float32_t b
	np.float32_t phi_rot
	
@cython.boundscheck(False)
@cython.wraparound(False)
def build_volume_phantom_2d_params(int phantom_size, object_2d[:] obj_params):
	"""
	build_volume_phantom_2d (model_parameters_filename,model_id, phantom_size)
	
	Takes in a input model_id and phantom_size and returns a phantom of phantom_size x phantom_size of type float32 numpy array.
	
	param: model_parameters_filename -- filename for the model parameters
	param: model_id -- a model id from the functions file
	param: phantom_size -- a phantom size in each dimension.
	
	returns: numpy float32 phantom array
	
	"""
	cdef Py_ssize_t i
	cdef np.ndarray[np.float32_t, ndim=2, mode="c"] phantom = np.zeros([phantom_size, phantom_size], dtype='float32')
	cdef float ret_val
	for i in range(obj_params.shape[0]):
		ret_val = buildPhantom2D_core_single(&phantom[0,0], phantom_size, obj_params[i].Obj, obj_params[i].C0, obj_params[i].x0, obj_params[i].y0, obj_params[i].a, obj_params[i].b, obj_params[i].phi_rot)
	return phantom	
	
@cython.boundscheck(False)
@cython.wraparound(False)
def buildPhantom2D(int model_id, int phantom_size, str model_parameters_filename):
	"""
	buildPhantom2D(model_id, phantom_size,model_parameters_filename)
	
	Takes in a input model_id and phantom_size and returns a phantom of phantom_size x phantom_size of type float32 numpy array.
	
	param: model_parameters_filename -- filename for the model parameters
	param: model_id -- a model id from the functions file
	param: phantom_size -- a phantom size in each dimension.
	
	returns: numpy float32 phantom array
	
	"""
	
	cdef np.ndarray[np.float32_t, ndim=2, mode="c"] phantom = np.zeros([phantom_size, phantom_size], dtype='float32')
	cdef float ret_val
	py_byte_string = model_parameters_filename.encode('UTF-8')
	cdef char* c_string = py_byte_string
	ret_val = buildPhantom2D_core(&phantom[0,0], model_id, phantom_size, c_string)
	return phantom

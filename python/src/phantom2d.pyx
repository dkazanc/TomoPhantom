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
cdef extern float buildSino2D_core(float *A, int ModelSelected, int N, int P, float *Th, int AngTot, int CenTypeIn, char* ModelParametersFilename)
cdef extern float buildSino2D_core_single(float *A, int N, int P, float *Th, int AngTot, int CenTypeIn, int Object, float C0, float x0, float y0, float a, float b, float phi_rot)
	
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

@cython.boundscheck(False)
@cython.wraparound(False)
def buildSino2D(int model_id, int image_size, int detector_size, np.ndarray[np.float32_t, ndim=1, mode="c"] angles, str model_parameters_filename, int CenTypeIn):
    
	"""
	buildSino2D (model_id, image_size, detector_size, angles, model_parameters_filename, CenTypeIn)
	
	Takes in as input model_id, image_size, detector_size and projection angles and return a 3D sinogram corresponding to the model id.
	
	param: model_parameters_filename -- filename for the model parameters
	param: model_id -- a model id from the functions file
	param: image_size -- a phantom size in each dimension.
	param: detector_size -- int detector size.
	param: angles -- a numpy array of float values with angles in radians
	param: CenTypeIn -- 1 as default [0: radon, 1:astra]
	returns: numpy float32 phantom sinograms array.
	
	"""
	cdef np.ndarray[np.float32_t, ndim=2, mode="c"] sinogram = np.zeros([angles.shape[0], detector_size], dtype='float32')
	cdef float ret_val
	py_byte_string = model_parameters_filename.encode('UTF-8')
	cdef char* c_string = py_byte_string    
	cdef int AngTot = angles.shape[0]
	ret_val = buildSino2D_core(&sinogram[0,0], model_id, image_size, detector_size, &angles[0], AngTot, CenTypeIn, c_string)
	return sinogram	
	
@cython.boundscheck(False)
@cython.wraparound(False)
def build_sinogram_phantom_2d_params(int image_size, int detector_size, np.ndarray[np.float32_t, ndim=1, mode="c"] angles, int CenTypeIn, object_2d[:] obj_params):
    
	"""
	build_sinogram_phantom_2d  (model_id, image_size, detector_size, angles, model_parameters_filename, CenTypeIn)
	
	Takes in as input model parameters list, image_size, detector_size and projection angles and return a 2D sinogram corresponding to the model id.
	
	param: image_size -- a phantom size in each dimension.
	param: detector_size -- int detector size.
	param: angles -- a numpy array of float values with angles in radians
	param: CenTypeIn -- 1 as default [0: radon, 1:astra]
	param: obj_params -- object parameters list
	returns: numpy float32 phantom sinograms array.
	
	"""
	cdef Py_ssize_t i	
	cdef np.ndarray[np.float32_t, ndim=2, mode="c"] sinogram = np.zeros([angles.shape[0], detector_size], dtype='float32')
	cdef float ret_val 
	cdef int AngTot = angles.shape[0]
	for i in range(obj_params.shape[0]):
		ret_val = buildSino2D_core_single(&sinogram[0,0], image_size, detector_size, &angles[0], AngTot, CenTypeIn, obj_params[i].Obj, obj_params[i].C0, obj_params[i].x0, obj_params[i].y0, obj_params[i].a, obj_params[i].b, obj_params[i].phi_rot)
	return sinogram		
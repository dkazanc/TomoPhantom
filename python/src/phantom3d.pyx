"""
tomophantom.pyx
Phantom 3D Generator
"""

import cython

# import numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern float buildPhantom3D_core(float *A, int ModelSelected, int N, char* ModelParametersFilename)

@cython.boundscheck(False)
@cython.wraparound(False)
def build_phantom_3d(str model_parameters_filename, int model_id, int phantom_size):
	"""
	build_phantom_3d (model_id, phantom_size)
	
	Takes in a input model_id and phantom_size and returns a phantom of phantom_size x phantom_size x phantom_size of type float32 numpy array.
	
	param: model_parameters_filename -- filename for the model parameters
	param: model_id -- a model id from the functions file
	param: phantom_size -- a phantom size in each dimension.
	
	returns: numpy float32 phantom array
	
	"""
	
	cdef np.ndarray[np.float32_t, ndim=3, mode="c"] phantom = np.empty([phantom_size, phantom_size, phantom_size], dtype='float32')
	cdef float ret_val
	py_byte_string = model_parameters_filename.encode('UTF-8')
	cdef char* c_string = py_byte_string
	ret_val = buildPhantom3D_core(&phantom[0,0,0], model_id, phantom_size, c_string)
	return phantom
	
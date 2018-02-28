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

import cython

# import numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# declare the interface to the C code
cdef extern float TomoP3DModel_core(float *A, int ModelSelected, int N, char* ModelParametersFilename)
cdef extern float TomoP3DObject(float *A, int N, char *Object, float C0, float x0, float y0, float z0, float a, float b, float c, float psi1, float psi2, float psi3)
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
	
@cython.boundscheck(False)
@cython.wraparound(False)
def Object(int phantom_size, object_3d[:] obj_params):
	"""
	Object (model_parameters_filename,model_id, phantom_size)
	
	Takes in a input object description (list) and phantom_size and returns a phantom-object (3D) of phantom_size x phantom_size x phantom_size of type float32 numpy array.
	
	param: obj_params -- object parameters list
	param: phantom_size -- a phantom size in each dimension.
	
	returns: numpy float32 phantom array
	
	"""
	cdef Py_ssize_t i
	cdef np.ndarray[np.float32_t, ndim=3, mode="c"] phantom = np.zeros([phantom_size, phantom_size, phantom_size], dtype='float32')
	cdef float ret_val
	for i in range(obj_params.shape[0]):
		ret_val = TomoP3DObject(&phantom[0,0,0], phantom_size, obj_params[i].Obj, obj_params[i].C0, obj_params[i].x0, obj_params[i].y0, obj_params[i].z0, obj_params[i].a, obj_params[i].b, obj_params[i].c, obj_params[i].psi1, obj_params[i].psi2, obj_params[i].psi3)
	return phantom	
	
@cython.boundscheck(False)
@cython.wraparound(False)
def Model(int model_id, int phantom_size, str model_parameters_filename):
	"""
	Model(model_id, phantom_size,model_parameters_filename)
	
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
	ret_val = TomoP3DModel_core(&phantom[0,0,0], model_id, phantom_size, c_string)
	return phantom

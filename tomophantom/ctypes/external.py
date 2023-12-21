#!/usr/bin/env python
# -*- coding: utf-8 -*-

#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
Ctypes-wrapped C-modules of TomoPhantom library.
"""

import tomophantom.ctypes.dtype as dtype
from . import c_shared_lib

__all__ = [
    "c_checkParams2D",
    "c_model2d",
    "c_model_sino2d",
    "c_object_sino2d",
    "c_object2d",
    "c_checkParams3D",
    "c_model3d",
    "c_object3d",
    "c_model_sino3d",
    "c_object_sino3d",
]

LIB_TOMOPHANTOM = c_shared_lib("tomophantom")


def c_checkParams2D(
    params,
    model_no,
    models_library_path,
):
    LIB_TOMOPHANTOM.checkParams2D.restype = dtype.as_c_void_p()
    LIB_TOMOPHANTOM.checkParams2D(
        dtype.as_c_int_p(params),
        dtype.as_c_int(model_no),
        dtype.as_c_char_p(models_library_path),
    )
    return params


def c_checkParams3D(
    params,
    model_no,
    models_library_path,
):
    LIB_TOMOPHANTOM.checkParams3D.restype = dtype.as_c_void_p()
    LIB_TOMOPHANTOM.checkParams3D(
        dtype.as_c_int_p(params),
        dtype.as_c_int(model_no),
        dtype.as_c_char_p(models_library_path),
    )
    return params


def c_numerical_sino2d(
    sinogram,
    input,
    phantom_size,
    detector_X,
    angles,
    angles_total,
):
    LIB_TOMOPHANTOM.TomoP2DSinoNum_core.restype = dtype.as_c_void_p()
    LIB_TOMOPHANTOM.TomoP2DSinoNum_core(
        dtype.as_c_float_p(sinogram),
        dtype.as_c_float_p(input),
        dtype.as_c_int(phantom_size),
        dtype.as_c_int(detector_X),
        dtype.as_c_float_p(angles),
        dtype.as_c_int(angles_total),
        0,
    )
    return sinogram


def c_model2d(
    output,
    model_no,
    phantom_size,
    models_library_path,
):
    LIB_TOMOPHANTOM.TomoP2DModel_core.restype = dtype.as_c_void_p()
    LIB_TOMOPHANTOM.TomoP2DModel_core(
        dtype.as_c_float_p(output),
        dtype.as_c_int(model_no),
        dtype.as_c_int(phantom_size),
        dtype.as_c_char_p(models_library_path),
    )
    return output


def c_model_sino2d(
    sinogram,
    model_no,
    phantom_size,
    detector_size,
    angles,
    angles_total,
    models_library_path,
):
    LIB_TOMOPHANTOM.TomoP2DModelSino_core.restype = dtype.as_c_void_p()
    LIB_TOMOPHANTOM.TomoP2DModelSino_core(
        dtype.as_c_float_p(sinogram),
        dtype.as_c_int(model_no),
        dtype.as_c_int(phantom_size),
        dtype.as_c_int(detector_size),
        dtype.as_c_float_p(angles),
        dtype.as_c_int(angles_total),
        1,
        dtype.as_c_char_p(models_library_path),
    )
    return sinogram


def c_object_sino2d(
    sinogram,
    phantom_size,
    detector_size,
    angles,
    angles_total,
    objectName,
    C0,
    y0,
    x0,
    a,
    b,
    phi_ang,
    tt,
):
    LIB_TOMOPHANTOM.TomoP2DObjectSino_core.restype = dtype.as_c_void_p()
    LIB_TOMOPHANTOM.TomoP2DObjectSino_core(
        dtype.as_c_float_p(sinogram),
        dtype.as_c_int(phantom_size),
        dtype.as_c_int(detector_size),
        dtype.as_c_float_p(angles),
        dtype.as_c_int(angles_total),
        1,
        dtype.as_c_char_p(objectName),
        dtype.as_c_float(C0),
        dtype.as_c_float(y0),
        dtype.as_c_float(x0),
        dtype.as_c_float(a),
        dtype.as_c_float(b),
        dtype.as_c_float(phi_ang),
        dtype.as_c_int(tt),
    )
    return sinogram


def c_object2d(object2d, phantom_size, objectName, C0, y0, x0, a, b, phi_ang, tt):
    LIB_TOMOPHANTOM.TomoP2DObject_core.restype = dtype.as_c_void_p()
    LIB_TOMOPHANTOM.TomoP2DObject_core(
        dtype.as_c_float_p(object2d),
        dtype.as_c_int(phantom_size),
        dtype.as_c_char_p(objectName),
        dtype.as_c_float(C0),
        dtype.as_c_float(x0),
        dtype.as_c_float(y0),
        dtype.as_c_float(b),
        dtype.as_c_float(a),
        dtype.as_c_float(phi_ang),
        dtype.as_c_int(tt),
    )
    return object2d


def c_model3d(
    output,
    model_no,
    N1,
    N2,
    N3,
    Z1,
    Z2,
    models_library_path,
):
    LIB_TOMOPHANTOM.TomoP3DModel_core.restype = dtype.as_c_void_p()
    LIB_TOMOPHANTOM.TomoP3DModel_core(
        dtype.as_c_float_p(output),
        dtype.as_c_int(model_no),
        dtype.as_c_long(N3),
        dtype.as_c_long(N2),
        dtype.as_c_long(N1),
        dtype.as_c_long(Z1),
        dtype.as_c_long(Z2),
        dtype.as_c_char_p(models_library_path),
    )
    return output

def c_object3d(
    output,
    N1,
    N2,
    N3,
    Z1,
    Z2,
    objectName,
    C0,
    x0,
    y0,
    z0,
    a,
    b,
    c,
    phi1,
    phi2,
    phi3,
    tt,
):
    LIB_TOMOPHANTOM.TomoP3DObject_core.restype = dtype.as_c_void_p()
    LIB_TOMOPHANTOM.TomoP3DObject_core(
        dtype.as_c_float_p(output),
        dtype.as_c_long(N3),
        dtype.as_c_long(N2),
        dtype.as_c_long(N1),
        dtype.as_c_long(Z1),
        dtype.as_c_long(Z2),
        dtype.as_c_char_p(objectName),
        dtype.as_c_float(C0),
        dtype.as_c_float(y0),
        dtype.as_c_float(x0),
        dtype.as_c_float(z0),
        dtype.as_c_float(a),
        dtype.as_c_float(b),
        dtype.as_c_float(c),
        dtype.as_c_float(phi1),
        dtype.as_c_float(phi2),
        dtype.as_c_float(phi3),
        dtype.as_c_long(tt),
    )
    return output

def c_model_sino3d(
    output,
    model_no,
    detector_horiz,
    detector_vert,
    Z1,
    Z2,
    phantom_size,
    angles,
    angles_total,
    models_library_path,
):
    LIB_TOMOPHANTOM.TomoP3DModelSino_core.restype = dtype.as_c_void_p()
    LIB_TOMOPHANTOM.TomoP3DModelSino_core(
        dtype.as_c_float_p(output),
        dtype.as_c_int(model_no),
        dtype.as_c_long(detector_horiz),
        dtype.as_c_long(detector_vert),        
        dtype.as_c_long(Z1),
        dtype.as_c_long(Z2),
        dtype.as_c_long(phantom_size),
        dtype.as_c_float_p(angles),
        dtype.as_c_int(angles_total),
        dtype.as_c_char_p(models_library_path),
    )
    return output

def c_object_sino3d(
    output,
    detector_horiz,
    detector_vert,
    Z1,
    Z2,
    phantom_size,
    angles,
    angles_total,
    objectName,
    C0,
    x0,
    y0,
    z0,
    a,
    b,
    c,
    phi1,
    phi2,
    phi3,
    tt,
):
    LIB_TOMOPHANTOM.TomoP3DObjectSino_core.restype = dtype.as_c_void_p()
    if (("gaussian" in str(objectName)) or ("paraboloid" in str(objectName)) or ("ellipsoid" in str(objectName))):
        LIB_TOMOPHANTOM.TomoP3DObjectSino_core(
            dtype.as_c_float_p(output),
            dtype.as_c_long(detector_horiz),
            dtype.as_c_long(detector_vert),
            dtype.as_c_long(Z1),
            dtype.as_c_long(Z2),
            dtype.as_c_long(phantom_size),
            dtype.as_c_float_p(angles),
            dtype.as_c_int(angles_total),
            dtype.as_c_char_p(objectName),
            dtype.as_c_float(C0),
            dtype.as_c_float(y0),
            dtype.as_c_float(-z0),
            dtype.as_c_float(-x0),
            dtype.as_c_float(b),
            dtype.as_c_float(a),
            dtype.as_c_float(c),
            dtype.as_c_float(phi3),
            dtype.as_c_float(phi2),
            dtype.as_c_float(phi1),
            dtype.as_c_long(tt),
        )
    elif ("elliptical_cylinder" in str(objectName)):
        LIB_TOMOPHANTOM.TomoP3DObjectSino_core(
            dtype.as_c_float_p(output),
            dtype.as_c_long(detector_horiz),
            dtype.as_c_long(detector_vert),
            dtype.as_c_long(Z1),
            dtype.as_c_long(Z2),
            dtype.as_c_long(phantom_size),
            dtype.as_c_float_p(angles),
            dtype.as_c_int(angles_total),
            dtype.as_c_char_p(objectName),
            dtype.as_c_float(C0),
            dtype.as_c_float(x0),
            dtype.as_c_float(-y0),
            dtype.as_c_float(z0),
            dtype.as_c_float(b),
            dtype.as_c_float(a),
            dtype.as_c_float(c),
            dtype.as_c_float(phi3),
            dtype.as_c_float(phi2),
            dtype.as_c_float(phi1),
            dtype.as_c_long(tt),
        )
    else:
        LIB_TOMOPHANTOM.TomoP3DObjectSino_core(
            dtype.as_c_float_p(output),
            dtype.as_c_long(detector_horiz),
            dtype.as_c_long(detector_vert),
            dtype.as_c_long(Z1),
            dtype.as_c_long(Z2),
            dtype.as_c_long(phantom_size),
            dtype.as_c_float_p(angles),
            dtype.as_c_int(angles_total),
            dtype.as_c_char_p(objectName),
            dtype.as_c_float(C0),
            dtype.as_c_float(x0),
            dtype.as_c_float(y0),
            dtype.as_c_float(z0),
            dtype.as_c_float(a),
            dtype.as_c_float(b),
            dtype.as_c_float(c),
            dtype.as_c_float(phi3),
            dtype.as_c_float(phi2),
            dtype.as_c_float(-phi1),
            dtype.as_c_long(tt),
        )
    return output
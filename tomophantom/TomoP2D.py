"""
Copyright 2023
The University of Manchester & Diamond Light Source
Licensed under the Apache License, Version 2.0.


These are modules for generation of 2D and 3D (dynamic extensions)
phantoms and their analytical data.

API Summary:

* Model - generates a 2D phantom from the Library (Phantom2DLibrary)
* ModelTemporal - generates a (2D+t) 3D model from the Library (Phantom2DLibrary)
* ModelSino - generates 2D sinogram for a model from the Library (Phantom2DLibrary)
* ModelSinoTemporal - generates (2D+t) sinogram for a temporal model from the Library (Phantom2DLibrary)
* Object - generates a 2D phantom from a given object described in the dictionary.
* ObjectSino - generates a 2D sinogram from a given object described in the dictionary.
* SinoNum - calculates a sinogram numerically from any given 2D array.
"""

import ctypes
import numpy as np
from numbers import Number
from enum import Enum
from typing import Union
from pathlib import Path

import tomophantom.ctypes.external as external

__all__ = [
    "Model",
    "ModelTemporal",
    "ModelSino",
    "ModelSinoTemporal",
    "Object",
    "ObjectSino",
    "SinoNum",
]


class Objects2D(Enum):
    """Enumeration with the available objects for 2D phantoms"""

    GAUSSIAN = "gaussian"
    PARABOLA = "parabola"
    PARABOLA1 = "parabola1"
    ELLIPSE = "ellipse"
    CONE = "cone"
    RECTANGLE = "rectangle"


def _check_params2d(model_no: int, models_library_path: Path) -> np.ndarray:
    """Check parameters before executing the generation script.

    Args:
        model_no (int): Model number from the Phantom2DLibrary.dat library file.
        models_library_path (Path): A path to the library file.

    Returns:
        list: a list of integers
    """
    params = np.ascontiguousarray(np.zeros([10], dtype=ctypes.c_int))
    external.c_checkParams2D(
        params,
        model_no,
        models_library_path,
    )
    return params


def Model(model_no: int, phantom_size: int, models_library_path: Path) -> np.ndarray:
    """Generate 2D phantoms based on the model no. in
    the library file Phantom2DLibrary.dat.

    Args:
        model_no (int): Model number from the Phantom2DLibrary.dat library file.
        phantom_size (int): A size of the generated phantom (squared).
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 2D phantom (N x N).
    """
    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([10], dtype=ctypes.c_int))
    params = _check_params2d(model_no, models_library_path)
    __testParams2D(params)  # check parameters and terminate if incorrect

    if params[3] == 1:
        # execute the model building function
        phantom = np.zeros([phantom_size, phantom_size], dtype="float32", order="C")
    else:
        raise ValueError(
            "The selected model is temporal (2D+time), use 'ModelTemporal' function instead"
        )
    external.c_model2d(
        np.ascontiguousarray(phantom),
        model_no,
        phantom_size,
        models_library_path,
    )
    return phantom


def ModelTemporal(
    model_no: int, phantom_size: int, models_library_path: Path
) -> np.ndarray:
    """Generate 2D+time temporal phantoms based on the model no. in
    the library file Phantom2DLibrary.dat. Note that temporal phantom
    numbers begin at 100 onwards.

    Args:
        model_no (int): Temporal model number from the Phantom2DLibrary.dat library file.
        phantom_size (int): A size of the generated phantom (squared).
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 2D+time phantom (time_frames x N x N).
    """
    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([10], dtype=ctypes.c_int))
    params = _check_params2d(model_no, models_library_path)
    __testParams2D(params)  # check parameters and terminate if incorrect

    if params[3] == 1:
        raise ValueError(
            "The selected model is stationary (2D), use 'Model' function instead"
        )
    else:
        # execute the model building function
        phantom = np.zeros(
            [params[3], phantom_size, phantom_size], dtype="float32", order="C"
        )
    external.c_model2d(
        np.ascontiguousarray(phantom),
        model_no,
        phantom_size,
        models_library_path,
    )
    return phantom


def ModelSino(
    model_no: int,
    phantom_size: int,
    detector_size: int,
    angles: np.ndarray,
    models_library_path: Path,
) -> np.ndarray:
    """Generate 2D analytical sinogram for corresponding models in
    the library file Phantom2DLibrary.dat.

    Args:
        model_no (int): Model number from the Phantom2DLibrary.dat library file.
        phantom_size (int): A size of the phantom (squared).
        detector_size (int): A size of the horizontal detector.
        angles (np.ndarray): Angles vector in degrees.
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 2D analytical sinogram (angles_total x detector_size).
    """
    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([10], dtype=ctypes.c_int))
    params = _check_params2d(model_no, models_library_path)
    __testParams2D(params)  # check parameters and terminate if incorrect

    angles_total = len(angles)
    if params[3] == 1:
        # execute the model building function
        sino2d = np.zeros([detector_size, angles_total], dtype="float32", order="C")
    else:
        raise ValueError(
            "The selected model is temporal (2D+time), use 'ModelSinoTemporal' function instead"
        )
    external.c_model_sino2d(
        np.ascontiguousarray(sino2d),
        model_no,
        phantom_size,
        detector_size,
        np.ascontiguousarray(angles),
        angles_total,
        models_library_path,
    )
    return sino2d.transpose()


def ModelSinoTemporal(
    model_no: int,
    phantom_size: int,
    detector_size: int,
    angles: np.ndarray,
    models_library_path: Path,
) -> np.ndarray:
    """Generate 2D+time (temporal )analytical sinogram for
    corresponding models in the library file Phantom2DLibrary.dat.

    Args:
        model_no (int): Temporal model number from the Phantom2DLibrary.dat library file.
        phantom_size (int): A size of the phantom (squared).
        detector_size (int): A size of the horizontal detector.
        angles (np.ndarray): Angles vector in degrees.
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 2D+time analytical sinogram (time_frames x detector_size x angles_total).
    """
    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([10], dtype=ctypes.c_int))
    params = _check_params2d(model_no, models_library_path)
    __testParams2D(params)  # check parameters and terminate if incorrect

    angles_total = len(angles)
    if params[3] == 1:
        raise ValueError(
            "The selected model is stationary (2D), use 'ModelSino' function instead"
        )
    else:
        # execute the model building function
        sino2d_t = np.zeros(
            [params[3], detector_size, angles_total], dtype="float32", order="C"
        )
    external.c_model_sino2d(
        np.ascontiguousarray(sino2d_t),
        model_no,
        phantom_size,
        detector_size,
        np.ascontiguousarray(angles),
        angles_total,
        models_library_path,
    )
    return sino2d_t


def Object(phantom_size: int, obj_params: Union[list, dict]) -> np.ndarray:
    """Generate a 2D analytical phantom for the standalone
       geometrical object that is parametrised in the "obj_params" dictionary.

    Args:
        phantom_size (int): A size of the phantom (squared).
        obj_params (list of dict or a dict): A dictionary with parameters of an object, see Demos.

    Returns:
        np.ndarray: The generated 2D analytical object.
    """
    if type(obj_params) is dict:
        obj_params = [obj_params]

    object2d = np.zeros([phantom_size, phantom_size], dtype="float32", order="C")
    # unpacking obj_params dictionary
    for obj in obj_params:
        if __testParams(obj):
            objectName = obj["Obj"].value
            C0 = obj["C0"]
            x0 = obj["x0"]
            y0 = obj["y0"]
            a = obj["a"]
            b = obj["b"]
            phi = obj["phi"]
            external.c_object2d(
                np.ascontiguousarray(object2d),
                phantom_size,
                objectName,
                C0,
                x0,
                y0,
                a,
                b,
                phi_ang=phi,
                tt=0,
            )
    return object2d


def ObjectSino(
    phantom_size: int,
    detector_size: int,
    angles: np.ndarray,
    obj_params: Union[list, dict],
) -> np.ndarray:
    """Generate a 2D analytical sinogram for the standalone
       geometrical object that is parametrised in the "obj_params" dictionary.

    Args:
        phantom_size (int): A size of the phantom (squared).
        detector_size (int): A size of the horizontal detector.
        angles (np.ndarray): Angles vector in degrees.
        obj_params (list of dicts or a dict): A dictionary with parameters of an object, see Demos.

    Returns:
        np.ndarray: The generated 2D analytical sinogram of an object.
    """
    angles_total = len(angles)
    sino2d = np.zeros([detector_size, angles_total], dtype="float32", order="C")
    # unpacking obj_params dictionary
    for obj in obj_params:
        if __testParams(obj):
            objectName = obj["Obj"].value
            C0 = obj["C0"]
            x0 = obj["x0"]
            y0 = obj["y0"]
            a = obj["a"]
            b = obj["b"]
            phi = obj["phi"]
            external.c_object_sino2d(
                np.ascontiguousarray(sino2d),
                phantom_size,
                detector_size,
                np.ascontiguousarray(angles),
                angles_total,
                objectName,
                C0,
                x0,
                y0,
                a,
                b,
                phi_ang=-phi,
                tt=0,
            )
    return sino2d.transpose()


def SinoNum(input: np.ndarray, detector_X: int, angles: np.ndarray) -> np.ndarray:
    """Numerical calculation of 2D sinogram from the 2D input.

    Args:
        input (np.ndarray): 2D object (e.g. a phantom).
        detector_X (int): The size of the detector X.
        angles (np.ndarray): Angles vector in degrees.

    Raises:
        ValueError: if input is 3D

    Returns:
        np.ndarray: Numerical sinogram calculated from the object
    """
    if np.ndim(input) == 3:
        raise ValueError("The accepted inputs must be 2D")
    phantom_size = input.shape[0]
    angles_total = len(angles)

    sinogram = np.zeros([detector_X, angles_total], dtype="float32", order="C")

    external.c_numerical_sino2d(
        np.ascontiguousarray(sinogram),
        np.ascontiguousarray(input),
        phantom_size,
        detector_X,
        np.ascontiguousarray(angles),
        angles_total,
    )
    return sinogram.transpose()


def __testParams(obj):
    """Performs a simple type check of the input parameters and a range check"""
    if not type(obj) is dict:
        raise TypeError("obj is not a dict {0}".format(type(obj)))

    # type check
    for k, v in obj.items():
        if not isinstance(v, Number):
            if not k == "Obj":
                raise TypeError(k, "is not a Number")
    typecheck = True
    # range check
    rangecheck = obj["x0"] >= -1 and obj["x0"] <= 1
    if not rangecheck:
        raise ValueError("x0 is out of range. Must be between -1 and 1")
    rangecheck = rangecheck and obj["y0"] >= -1 and obj["y0"] <= 1
    if not rangecheck:
        raise ValueError("y0 is out of range. Must be between -1 and 1")
    rangecheck = rangecheck and obj["a"] > 0 and obj["a"] <= 2
    if not rangecheck:
        raise ValueError("a (object size) must be positive in [0,2] range")
    rangecheck = rangecheck and obj["b"] > 0 and obj["b"] <= 2
    if not rangecheck:
        raise ValueError("b (object size) must be positive in [0,2] range")
    return rangecheck and typecheck


def __testParams2D(obj):
    if obj[0] == 0:
        raise TypeError(
            "Check if the library file <Phantom2DLibrary.dat> exists and the given path is correct"
        )
    if obj[1] == 0:
        raise TypeError(
            "The given model is not found, check available models in <Phantom2DLibrary.dat> file"
        )
    if obj[2] == 0:
        raise TypeError(
            "Components number cannot be negative, check <Phantom2DLibrary.dat> file"
        )
    if obj[3] == 0:
        raise TypeError(
            "TimeSteps cannot be negative, check <Phantom2DLibrary.dat> file"
        )
    if obj[4] == 0:
        raise TypeError("Unknown name of the object, check <Phantom2DLibrary.dat> file")
    if obj[5] == 0:
        raise TypeError(
            "C0 should not be equal to zero, check <Phantom2DLibrary.dat> file"
        )
    if obj[6] == 0:
        raise TypeError(
            "x0 (object position) must be in [-1,1] range, check <Phantom2DLibrary.dat> file"
        )
    if obj[7] == 0:
        raise TypeError(
            "y0 (object position) must be in [-1,1] range, check <Phantom2DLibrary.dat> file"
        )
    if obj[8] == 0:
        raise TypeError(
            "a (object size) must be positive in [0,2] range, check <Phantom2DLibrary.dat> file"
        )
    if obj[9] == 0:
        raise TypeError(
            "b (object size) must be positive in [0,2] range, check <Phantom2DLibrary.dat> file"
        )
    return 0

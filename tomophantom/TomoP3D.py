"""
Copyright 2023
The University of Manchester & Diamond Light Source
Licensed under the Apache License, Version 2.0.


These are modules for generation of 3D and 4D (dynamic extensions) phantoms and
their analytical data.

API Summary:

* Model - generates a 3D phantom from the Library (Phantom3DLibrary)
* ModelSub - generates a vertical subset (cutoff) for model from the Library (Phantom3DLibrary)
* ModelTemporal - generates a (3D + t) 4D model from the Library (Phantom3DLibrary)
* ModelTemporalSub - generates a vertical subset (cutoff) of 4D model from the Library (Phantom3DLibrary)
* ModelSino - generates 3D projection data for a model from the Library (Phantom3DLibrary)
* ModelSinoSub - generates a vertical subset (cutoff) for 3D projection data for a model from the Library (Phantom3DLibrary)
* ModelSinoTemporal - generates (3D + t) 4D projection data for a temporal model from the Library (Phantom3DLibrary)
* ModelSinoTemporalSub - generates (3D + t) 4D projection data (subset) for a temporal model from the Library (Phantom3DLibrary)
* Object - generates a 3D phantom from a given object described in the dictionary.
* ObjectSino - generates a 3D projection data from a given object described in the dictionary.
"""

import ctypes
import numpy as np
from numbers import Number
from enum import Enum
from typing import Union, Tuple
from pathlib import Path

import tomophantom.ctypes.external as external

__all__ = [
    "Model",
    "ModelSub",
    "ModelTemporal",
    "ModelTemporalSub",
    "ModelSino",
    "ModelSinoSub",
    "ModelSinoTemporal",
    "ModelSinoTemporalSub",
    "Object",
    "ObjectSino",
]


class Objects3D(Enum):
    """Enumeration with the available objects for 3D phantoms"""

    GAUSSIAN = "gaussian"
    PARABOLOID = "paraboloid"
    ELLIPSOID = "ellipsoid"
    CONE = "cone"
    CUBOID = "cuboid"
    ELLIPCYLINDER = "elliptical_cylinder"


def _check_params3d(model_no: int, models_library_path: Path) -> np.ndarray:
    """Check parameters before executing the generation script.

    Args:
        model_no (int): Model number from the Phantom3DLibrary.dat library file.
        models_library_path (Path): A path to the library file.

    Returns:
        list: a list of integers
    """
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    external.c_checkParams3D(
        params,
        model_no,
        models_library_path,
    )
    return params


def Model(
    model_no: int,
    phantom_size: Union[int, Tuple[int, int, int]],
    models_library_path: Path,
) -> np.ndarray:
    """Generates a 3D phantom based on the model no. in the library file Phantom3DLibrary.dat.

    Args:
        model_no (int): Model number from the Phantom3DLibrary.dat library file.
        phantom_size (int, Tuple(int)): A scalar or a tuple with phantom dimensions.
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 3D phantom.
    """
    if type(phantom_size) == tuple:
        N1, N2, N3 = [int(i) for i in phantom_size]
    else:
        N1 = N2 = N3 = phantom_size

    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    params = _check_params3d(model_no, models_library_path)
    __testParams3D(params)  # check parameters and terminate if incorrect

    if params[3] == 1:
        phantom = np.zeros([N1, N2, N3], dtype="float32", order="C")
    else:
        raise ValueError(
            "The selected model is temporal (3D+time), use 'ModelTemporal' function instead"
        )
    external.c_model3d(
        np.ascontiguousarray(phantom),
        model_no,
        N1,
        N2,
        N3,
        0,
        N1,
        models_library_path,
    )
    return phantom


def ModelSub(
    model_no: int,
    phantom_size: Union[int, Tuple[int, int, int]],
    sub_index: Tuple[int, int],
    models_library_path: Path,
) -> np.ndarray:
    """Generates a subset of 3D phantom (vertical cutoff) based on the model no. in
    the library file Phantom3DLibrary.dat.

    Args:
        model_no (int): Model number from the Phantom3DLibrary.dat library file.
        phantom_size (int, Tuple(int)): A scalar or a tuple with phantom dimensions.
        sub_index (Tuple[int, int]): a tuple containing 2 indeces (lower, upper) to select a vertical subset.
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 3D phantom.
    """
    if type(phantom_size) == tuple:
        N1, N2, N3 = [int(i) for i in phantom_size]
    else:
        N1 = N2 = N3 = phantom_size

    Z1, Z2 = [int(i) for i in sub_index]

    rangecheck = Z2 > Z1
    if not rangecheck:
        raise ValueError("Upper index must be larger than the lower one")
    rangecheck = Z1 >= 0 and Z1 < N3
    if not rangecheck:
        raise ValueError("Range of the lower index is incorrect")
    rangecheck = Z2 >= 0 and Z2 <= N3
    if not rangecheck:
        raise ValueError("Range of the higher index is incorrect")

    sub_size = Z2 - Z1  # the size of the vertical slab

    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    params = _check_params3d(model_no, models_library_path)
    __testParams3D(params)  # check parameters and terminate if incorrect

    if params[3] == 1:
        phantom = np.zeros([sub_size, N2, N3], dtype="float32", order="C")
    else:
        raise ValueError(
            "The selected model is temporal (3D+time), use 'ModelTemporal' function instead"
        )
    external.c_model3d(
        np.ascontiguousarray(phantom),
        model_no,
        N1,
        N2,
        N3,
        Z1,
        Z2,
        models_library_path,
    )
    return phantom


def ModelTemporal(
    model_no: int,
    phantom_size: Union[int, Tuple[int, int, int]],
    models_library_path: Path,
) -> np.ndarray:
    """Generates 4D (3D+time) phantom based on the model no. in
    the library file Phantom3DLibrary.dat.

    Args:
        model_no (int): Model number from the Phantom3DLibrary.dat library file.
        phantom_size (int, Tuple(int)): A scalar or a tuple with phantom dimensions.
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 4D phantom.
    """
    if type(phantom_size) == tuple:
        N1, N2, N3 = [int(i) for i in phantom_size]
    else:
        N1 = N2 = N3 = phantom_size

    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    params = _check_params3d(model_no, models_library_path)
    __testParams3D(params)  # check parameters and terminate if incorrect

    if params[3] == 1:
        raise ValueError(
            "The selected model is stationary (3D), use 'Model' function instead"
        )
    else:
        phantom = np.zeros([params[3], N1, N2, N3], dtype="float32", order="C")
    external.c_model3d(
        np.ascontiguousarray(phantom),
        model_no,
        N1,
        N2,
        N3,
        0,
        N1,
        models_library_path,
    )
    return phantom


def ModelTemporalSub(
    model_no: int,
    phantom_size: Union[int, Tuple[int, int, int]],
    sub_index: Tuple[int, int],
    models_library_path: Path,
) -> np.ndarray:
    """Generates a subset of 4D phantom (vertical cutoff) based on the model no. in
    the library file Phantom3DLibrary.dat.

    Args:
        model_no (int): Model number from the Phantom3DLibrary.dat library file.
        phantom_size (int, Tuple(int)): A scalar or a tuple with phantom dimensions.
        sub_index (Tuple[int, int]): a tuple containing 2 indeces (lower, upper) to select a vertical subset.
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 4D phantom.
    """
    if type(phantom_size) == tuple:
        N1, N2, N3 = [int(i) for i in phantom_size]
    else:
        N1 = N2 = N3 = phantom_size

    Z1, Z2 = [int(i) for i in sub_index]

    rangecheck = Z2 > Z1
    if not rangecheck:
        raise ValueError("Upper index must be larger than the lower one")
    rangecheck = Z1 >= 0 and Z1 < N3
    if not rangecheck:
        raise ValueError("Range of the lower index is incorrect")
    rangecheck = Z2 >= 0 and Z2 <= N3
    if not rangecheck:
        raise ValueError("Range of the higher index is incorrect")

    sub_size = Z2 - Z1  # the size of the vertical slab

    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    params = _check_params3d(model_no, models_library_path)
    __testParams3D(params)  # check parameters and terminate if incorrect

    if params[3] == 1:
        raise ValueError(
            "The selected model is temporal (3D+time), use 'ModelTemporal' function instead"
        )
    else:
        phantom = np.zeros([params[3], sub_size, N2, N3], dtype="float32", order="C")
    external.c_model3d(
        np.ascontiguousarray(phantom),
        model_no,
        N1,
        N2,
        N3,
        Z1,
        Z2,
        models_library_path,
    )
    return phantom


def ModelSino(
    model_no: int,
    phantom_size: int,
    detector_horiz: int,
    detector_vert: int,
    angles: np.ndarray,
    models_library_path: Path,
) -> np.ndarray:
    """Generates a 3D analytical sinogram for corresponding models in
    the library file Phantom3DLibrary.dat.

    Args:
        model_no (int): Model number from the Phantom3DLibrary.dat library file.
        phantom_size (int): A scalar for squared phantom dimensions (2D slice).
        detector_horiz (int): Size of the horizontal detector in pixels.
        detector_vert (int): Size of the vertical detector in pixels.
        angles (np.ndarray): Angles vector in degrees.
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 3D projection data of the following shape [detector_vert, AngTot, detector_horiz].
    """
    if type(phantom_size) == tuple:
        raise ValueError(
            "Please give a scalar for the phantom size (2d slice), if one needs non-cubic data generator please see ModelSinoSub function"
        )

    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    params = _check_params3d(model_no, models_library_path)
    __testParams3D(params)  # check parameters and terminate if incorrect

    angles_total = len(angles)
    if params[3] == 1:
        # execute the model building function
        sino3d = np.zeros(
            [angles_total, detector_vert, detector_horiz], dtype="float32", order="C"
        )
    else:
        raise ValueError(
            "The selected model is temporal (4D), use 'ModelTemporalSino' function instead"
        )
    external.c_model_sino3d(
        np.ascontiguousarray(sino3d),
        model_no,
        detector_horiz,
        detector_vert,
        0,
        detector_vert,
        phantom_size,
        np.ascontiguousarray(angles),
        angles_total,
        models_library_path,
    )
    return np.swapaxes(sino3d, 0, 1)


def ModelSinoSub(
    model_no: int,
    phantom_size: int,
    detector_horiz: int,
    detector_vert: int,
    sub_index: tuple,
    angles: np.ndarray,
    models_library_path: Path,
) -> np.ndarray:
    """Generates a 3D analytical sinogram with vertical cutoff using sub_index tuple for corresponding models in
    the library file Phantom3DLibrary.dat.

    Args:
        model_no (int): Model number from the Phantom3DLibrary.dat library file.
        phantom_size (int): A scalar for squared phantom dimensions (2D slice).
        detector_horiz (int): Size of the horizontal detector in pixels.
        detector_vert (int): Size of the vertical detector in pixels.
        sub_index (Tuple[int, int]): a tuple containing 2 indeces (lower, upper) to select a vertical subset.
        angles (np.ndarray): Angles vector in degrees.
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 3D projection data of the following shape [detector_vert, AngTot, detector_horiz].
    """
    if type(phantom_size) == tuple:
        raise ValueError(
            "Please give a scalar for phantom size, projection data cannot be obtained for non-cubic phantom"
        )

    Z1, Z2 = [int(i) for i in sub_index]

    rangecheck = Z2 > Z1
    if not rangecheck:
        raise ValueError("Upper index must be larger than the lower one")
    rangecheck = Z1 >= 0 and Z1 < detector_vert
    if not rangecheck:
        raise ValueError("Range of the lower index is incorrect")
    rangecheck = Z2 >= 0 and Z2 <= detector_vert
    if not rangecheck:
        raise ValueError("Range of the higher index is incorrect")

    sub_size = Z2 - Z1  # the size of the vertical slab

    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    params = _check_params3d(model_no, models_library_path)
    __testParams3D(params)  # check parameters and terminate if incorrect

    angles_total = len(angles)
    if params[3] == 1:
        # execute the model building function
        sino3d = np.zeros(
            [angles_total, sub_size, detector_horiz], dtype="float32", order="C"
        )
    else:
        raise ValueError(
            "The selected model is temporal (4D), use 'ModelTemporalSino' function instead"
        )
    external.c_model_sino3d(
        np.ascontiguousarray(sino3d),
        model_no,
        detector_horiz,
        detector_vert,
        Z1,
        Z2,
        phantom_size,
        np.ascontiguousarray(angles),
        angles_total,
        models_library_path,
    )
    return np.swapaxes(sino3d, 0, 1)


def ModelSinoTemporal(
    model_no: int,
    phantom_size: int,
    detector_horiz: int,
    detector_vert: int,
    angles: np.ndarray,
    models_library_path: Path,
) -> np.ndarray:
    """Generates a 4D (3D+time) analytical sinogram for corresponding models in
    the library file Phantom3DLibrary.dat.

    Args:
        model_no (int): Temporal model number from the Phantom3DLibrary.dat library file.
        phantom_size (int): A scalar for squared phantom dimensions (2D slice).
        detector_horiz (int): Size of the horizontal detector in pixels.
        detector_vert (int): Size of the vertical detector in pixels.
        angles (np.ndarray): Angles vector in degrees.
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 4D projection data of the following shape [time_frames, detector_vert, AngTot, detector_horiz].
    """

    if type(phantom_size) == tuple:
        raise ValueError(
            "Please give a scalar for phantom size, projection data cannot be obtained for non-cubic phantom"
        )

    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    params = _check_params3d(model_no, models_library_path)
    __testParams3D(params)  # check parameters and terminate if incorrect

    angles_total = len(angles)
    if params[3] == 1:
        raise ValueError(
            "The selected model is stationary (3D), use 'ModelSino' function instead"
        )
    else:
        sino4d = np.zeros(
            [params[3], angles_total, detector_vert, detector_horiz],
            dtype="float32",
            order="C",
        )

    external.c_model_sino3d(
        np.ascontiguousarray(sino4d),
        model_no,
        detector_horiz,
        detector_vert,
        0,
        detector_vert,
        phantom_size,
        np.ascontiguousarray(angles),
        angles_total,
        models_library_path,
    )
    return np.swapaxes(sino4d, 1, 2)


def ModelSinoTemporalSub(
    model_no: int,
    phantom_size: int,
    detector_horiz: int,
    detector_vert: int,
    sub_index: tuple,
    angles: np.ndarray,
    models_library_path: Path,
) -> np.ndarray:
    """Generates a 4D analytical sinogram with vertical cutoff using sub_index tuple for corresponding models in
    the library file Phantom3DLibrary.dat.

    Args:
        model_no (int): Temporal model number from the Phantom3DLibrary.dat library file.
        phantom_size (int): A scalar for squared phantom dimensions (2D slice).
        detector_horiz (int): Size of the horizontal detector in pixels.
        detector_vert (int): Size of the vertical detector in pixels.
        sub_index (Tuple[int, int]): a tuple containing 2 indeces (lower, upper) to select a vertical subset.
        angles (np.ndarray): Angles vector in degrees.
        models_library_path (Path): A path to the library file.

    Returns:
        np.ndarray: The generated 3D projection data of the following shape [time, detector_vert, AngTot, detector_horiz].
    """
    if type(phantom_size) == tuple:
        raise ValueError(
            "Please give a scalar for phantom size, projection data cannot be obtained for non-cubic phantom"
        )

    Z1, Z2 = [int(i) for i in sub_index]

    rangecheck = Z2 > Z1
    if not rangecheck:
        raise ValueError("Upper index must be larger than the lower one")
    rangecheck = Z1 >= 0 and Z1 < detector_vert
    if not rangecheck:
        raise ValueError("Range of the lower index is incorrect")
    rangecheck = Z2 >= 0 and Z2 <= detector_vert
    if not rangecheck:
        raise ValueError("Range of the higher index is incorrect")

    sub_size = Z2 - Z1  # the size of the vertical slab

    # check the validity of model's parameters
    params = np.ascontiguousarray(np.zeros([12], dtype=ctypes.c_int))
    params = _check_params3d(model_no, models_library_path)
    __testParams3D(params)  # check parameters and terminate if incorrect

    angles_total = len(angles)
    if params[3] == 1:
        raise ValueError(
            "The selected model is stationary (3D), use 'ModelSino' function instead"
        )
    else:
        sino4d = np.zeros(
            [params[3], angles_total, sub_size, detector_horiz],
            dtype="float32",
            order="C",
        )

    external.c_model_sino3d(
        np.ascontiguousarray(sino4d),
        model_no,
        detector_horiz,
        detector_vert,
        Z1,
        Z2,
        phantom_size,
        np.ascontiguousarray(angles),
        angles_total,
        models_library_path,
    )
    return np.swapaxes(sino4d, 1, 2)


def Object(
    phantom_size: Union[int, Tuple[int, int, int]], obj_params: Union[list, dict]
) -> np.ndarray:
    """Generates a 3D object for the standalone geometrical figure
    that is parametrised in the "obj_params" dictionary.

    Args:
        phantom_size (int, Tuple(int)): A scalar or a tuple with phantom dimensions.
        obj_params (a list of dicts or dict): A dictionary with parameters of an object, see demos.

    Returns:
        np.ndarray: The generated 3D object.
    """

    if type(phantom_size) == tuple:
        N1, N2, N3 = [int(i) for i in phantom_size]
    else:
        N1 = N2 = N3 = phantom_size

    if type(obj_params) is dict:
        obj_params = [obj_params]

    object3d = np.zeros([N1, N2, N3], dtype="float32", order="C")
    # unpacking obj_params dictionary
    for obj in obj_params:
        if __testParams3d(obj):
            objectName = obj["Obj"].value
            C0 = obj["C0"]
            x0 = obj["x0"]
            y0 = obj["y0"]
            z0 = obj["z0"]
            a = obj["a"]
            b = obj["b"]
            c = obj["c"]
            phi1 = obj["phi1"]
            phi2 = 0.0
            phi3 = 0.0
            tt = 0
            external.c_object3d(
                np.ascontiguousarray(object3d),
                N1,
                N2,
                N3,
                0,
                N1,
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
            )
    return object3d


def ObjectSino(
    phantom_size: int,
    detector_horiz: int,
    detector_vert: int,
    angles: np.ndarray,
    obj_params: Union[list, dict],
) -> np.ndarray:
    """Generates a 3D analytical projection data for the standalone geometrical figure
    that is parametrised in the "obj_params" dictionary.

    Args:
        phantom_size (int): A scalar for phantom dimensions.
        detector_horiz (int): Size of the horizontal detector in pixels.
        detector_vert (int): Size of the vertical detector in pixels.
        angles (np.ndarray): Angles vector in degrees.
        obj_params (a list of dicts or dict): A dictionary with parameters of an object, see demos.

    Returns:
        np.ndarray: The generated 3D projection data for an object.
    """
    if type(phantom_size) == tuple:
        raise ValueError(
            "Please give a scalar for phantom size, projection data cannot be obtained for non-cubic phantom"
        )

    if type(obj_params) is dict:
        obj_params = [obj_params]

    angles_total = len(angles)
    tt = 0
    sino3d = np.zeros(
        [angles_total, detector_vert, detector_horiz], dtype="float32", order="C"
    )
    # unpacking obj_params dictionary
    for obj in obj_params:
        if __testParams3d(obj):
            objectName = obj["Obj"].value
            C0 = obj["C0"]
            x0 = obj["x0"]
            y0 = obj["y0"]
            z0 = obj["z0"]
            a = obj["a"]
            b = obj["b"]
            c = obj["c"]
            phi1 = obj["phi1"]
            phi2 = 0.0
            phi3 = 0.0
            tt = 0
            external.c_object_sino3d(
                np.ascontiguousarray(sino3d),
                detector_horiz,
                detector_vert,
                0,
                detector_vert,
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
            )
    return np.swapaxes(sino3d, 0, 1)


def __testParams3d(obj):
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
    rangecheck = rangecheck and obj["z0"] >= -1 and obj["z0"] <= 1
    if not rangecheck:
        raise ValueError("z0 is out of range. Must be between -1 and 1")
    rangecheck = rangecheck and obj["a"] > 0 and obj["a"] <= 2
    if not rangecheck:
        raise ValueError("a (object size) must be positive in [0,2] range")
    rangecheck = rangecheck and obj["b"] > 0 and obj["b"] <= 2
    if not rangecheck:
        raise ValueError("b (object size) must be positive in [0,2] range")
    rangecheck = rangecheck and obj["c"] > 0 and obj["c"] <= 2
    if not rangecheck:
        raise ValueError("c (object size) must be positive in [0,2] range")
    return rangecheck and typecheck


def __testParams3D(obj):
    if obj[0] == 0:
        raise TypeError(
            "Check if the library file <Phantom3DLibrary.dat> exists, the given path is correct and the syntax is valid"
        )
    if obj[1] == 0:
        raise TypeError(
            "The given model is not found, check available models in <Phantom3DLibrary.dat> file"
        )
    if obj[2] == 0:
        raise TypeError(
            "Components number cannot be negative, check <Phantom3DLibrary.dat> file"
        )
    if obj[3] == 0:
        raise TypeError(
            "TimeSteps cannot be negative, check <Phantom3DLibrary.dat> file"
        )
    if obj[4] == 0:
        raise TypeError("Unknown name of the object, check <Phantom3DLibrary.dat> file")
    if obj[5] == 0:
        raise TypeError(
            "C0 should not be equal to zero, check <Phantom3DLibrary.dat> file"
        )
    if obj[6] == 0:
        raise TypeError(
            "x0 (object position) must be in [-1,1] range, check <Phantom3DLibrary.dat> file"
        )
    if obj[7] == 0:
        raise TypeError(
            "y0 (object position) must be in [-1,1] range, check <Phantom3DLibrary.dat> file"
        )
    if obj[8] == 0:
        raise TypeError(
            "z0 (object position) must be in [-1,1] range, check <Phantom3DLibrary.dat> file"
        )
    if obj[9] == 0:
        raise TypeError(
            "a (object size) must be positive in [0,2] range, check <Phantom3DLibrary.dat> file"
        )
    if obj[10] == 0:
        raise TypeError(
            "b (object size) must be positive in [0,2] range, check <Phantom3DLibrary.dat> file"
        )
    if obj[11] == 0:
        raise TypeError(
            "c (object size) must be positive in [0,2] range, check <Phantom3DLibrary.dat> file"
        )
    return 0

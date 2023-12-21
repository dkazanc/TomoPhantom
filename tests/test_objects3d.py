import numpy as np
from numpy.testing import assert_allclose

from tomophantom.TomoP3D import _check_params3d, Objects3D, Object

eps = 1e-05


def test_3d_object():
    N_size = 64
    obj3D_1 = {
        "Obj": Objects3D.GAUSSIAN,
        "C0": 1.0,
        "x0": -0.25,
        "y0": -0.15,
        "z0": 0.0,
        "a": 0.3,
        "b": 0.2,
        "c": 0.3,
        "phi1": 35.0,
    }

    obj3D_2 = {
        "Obj": Objects3D.CUBOID,
        "C0": 1.00,
        "x0": 0.1,
        "y0": 0.2,
        "z0": 0.0,
        "a": 0.15,
        "b": 0.35,
        "c": 0.6,
        "phi1": -60.0,
    }

    myObjects = [obj3D_1, obj3D_2]  # a list of dictionaries
    Object3D = Object(N_size, myObjects)

    assert 0.001 <= np.max(Object3D) <= 100
    assert Object3D.dtype == np.float32
    assert Object3D.shape == (N_size, N_size, N_size)

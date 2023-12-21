import numpy as np
from numpy.testing import assert_allclose

import pytest
import os
import tomophantom
from tomophantom.TomoP2D import Objects2D, ObjectSino

eps = 1e-05


def test_2d_object_sino():
    N_size = 64  # set the desired dimension of the phantom
    angles_num = int(0.5 * np.pi * N_size)
    # angles number
    angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")
    P = int(np.sqrt(2) * N_size)  # detectors

    # define all objects bellow:
    pp = {
        "Obj": Objects2D.GAUSSIAN,
        "C0": 1.00,
        "x0": 0.25,
        "y0": -0.3,
        "a": 0.15,
        "b": 0.3,
        "phi": -30.0,
    }

    pp1 = {
        "Obj": Objects2D.RECTANGLE,
        "C0": 1.00,
        "x0": -0.2,
        "y0": 0.2,
        "a": 0.25,
        "b": 0.4,
        "phi": 60.0,
    }

    myObjects = [pp, pp1]  # list of dicts
    sino_an = ObjectSino(N_size, P, angles, myObjects)

    assert 0.0 <= np.max(sino_an) <= 1000
    assert sino_an.dtype == np.float32
    assert sino_an.shape == (angles_num, P)

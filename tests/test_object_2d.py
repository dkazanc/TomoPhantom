import numpy as np
from numpy.testing import assert_allclose

import os
import tomophantom
from tomophantom.TomoP2D import Objects2D, Object
from tomophantom.supp.libraryToDict import modelfile2Dtolist

eps = 1e-05


def test_2d_object():
    N_size = 64  # set the desired dimension of the phantom
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

    myObjects = [pp, pp1]  # dictionary of objects
    Object1 = Object(N_size, myObjects)

    assert 0.0 <= np.max(Object1) <= 3
    assert Object1.dtype == np.float32
    assert Object1.shape == (N_size, N_size)


def test_2d_object_extraction():
    model = 11  # select a model number from the library
    N_size = 64  # set dimension of the phantom
    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library2D = os.path.join(path, "phantomlib", "Phantom2DLibrary.dat")
    # extract parameters into a list
    objlist = modelfile2Dtolist(path_library2D, model)
    # This will generate a N_size x N_size phantom (2D)
    phantom_2D = Object(N_size, objlist)

    assert 0.0 <= np.max(phantom_2D) <= 3
    assert phantom_2D.dtype == np.float32
    assert phantom_2D.shape == (N_size, N_size)

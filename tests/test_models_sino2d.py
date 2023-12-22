import numpy as np
from numpy.testing import assert_allclose

import pytest
import os
import tomophantom
from tomophantom.TomoP2D import (
    _check_params2d,
    Model,
    ModelSino,
    SinoNum,
    ModelSinoTemporal,
)

eps = 1e-05


@pytest.mark.parametrize("model", list(np.arange(1, 16, 1)))
def test_2d_models_sino(model):
    N_size = 64  # set the desired dimension of the phantom
    angles_num = int(0.5 * np.pi * N_size)
    # angles number
    angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")
    P = int(np.sqrt(2) * N_size)  # detectors

    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library2D = os.path.join(path, "phantomlib", "Phantom2DLibrary.dat")

    # Generate a N_size x N_size phantom (2D)
    sino_an = ModelSino(model, N_size, P, angles, path_library2D)
    assert 0.0 <= np.max(sino_an) <= 1000
    assert sino_an.dtype == np.float32
    assert sino_an.shape == (angles_num, P)


def test_num_sino_2d():
    N_size = 64  # set the desired dimension of the phantom
    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library2D = os.path.join(path, "phantomlib", "Phantom2DLibrary.dat")

    angles_num = int(0.5 * np.pi * N_size)
    # angles number
    angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")
    P = int(np.sqrt(2) * N_size)  # detectors

    phantom_2D = Model(1, N_size, path_library2D)
    sino_num = SinoNum(phantom_2D, P, angles)

    assert sino_num.dtype == np.float32
    assert sino_num.shape == (100, 90)


@pytest.mark.parametrize("model", list(np.arange(100, 103, 1)))
def test_2d_temporal_sinos(model):
    N_size = 64  # set the desired dimension of the phantom
    angles_num = int(0.5 * np.pi * N_size)
    # angles number
    angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")
    P = int(np.sqrt(2) * N_size)  # detectors

    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library2D = os.path.join(path, "phantomlib", "Phantom2DLibrary.dat")

    params = _check_params2d(model, path_library2D)
    # Generate a time_frame x N_size x N_size phantom (3D)
    sino2d_t = ModelSinoTemporal(model, N_size, P, angles, path_library2D)
    assert 0.001 <= np.max(sino2d_t) <= 1000
    assert sino2d_t.dtype == np.float32
    assert sino2d_t.shape == (params[3], P, angles_num)

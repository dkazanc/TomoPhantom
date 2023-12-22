import numpy as np
from numpy.testing import assert_allclose

import pytest
import os
import tomophantom
from tomophantom.TomoP2D import _check_params2d, Model, ModelTemporal, SinoNum

eps = 1e-05


@pytest.mark.parametrize("model", list(np.arange(1, 16, 1)))
def test_2d_params_stat(model):
    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library2D = os.path.join(path, "phantomlib", "Phantom2DLibrary.dat")

    # check parameters
    params = _check_params2d(model, path_library2D)

    assert np.sum(params) == 10
    assert params.dtype == np.int32


@pytest.mark.parametrize("model", list(np.arange(100, 103, 1)))
def test_2d_params_temporal(model):
    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library2D = os.path.join(path, "phantomlib", "Phantom2DLibrary.dat")

    # check parameters
    params = _check_params2d(model, path_library2D)
    assert len(params) == 10
    if model == 100:
        assert params[3] == 3
    if model == 101:
        assert params[3] == 350
    if model == 102:
        assert params[3] == 25
    assert params.dtype == np.int32


@pytest.mark.parametrize("model", list(np.arange(1, 16, 1)))
def test_2d_models(model):
    N_size = 64  # set the desired dimension of the phantom
    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library2D = os.path.join(path, "phantomlib", "Phantom2DLibrary.dat")

    # Generate a N_size x N_size phantom (2D)
    phantom_2D = Model(model, N_size, path_library2D)
    assert 0.001 <= np.max(phantom_2D) <= 100
    assert phantom_2D.dtype == np.float32
    assert phantom_2D.shape == (N_size, N_size)


@pytest.mark.parametrize("model", list(np.arange(100, 103, 1)))
def test_2d_temporal_models(model):
    N_size = 64  # set the desired dimension of the phantom
    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library2D = os.path.join(path, "phantomlib", "Phantom2DLibrary.dat")

    params = _check_params2d(model, path_library2D)
    # Generate a time_frame x N_size x N_size phantom (3D)
    phantom_2D = ModelTemporal(model, N_size, path_library2D)
    assert 0.001 <= np.max(phantom_2D) <= 1000
    assert phantom_2D.dtype == np.float32
    assert phantom_2D.shape == (params[3], N_size, N_size)

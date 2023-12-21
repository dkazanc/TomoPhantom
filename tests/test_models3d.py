import numpy as np
from numpy.testing import assert_allclose

import pytest
import os
import tomophantom
from tomophantom.TomoP3D import (
    _check_params3d,
    Model,
    ModelSub,
    ModelTemporal,
    ModelTemporalSub,
)

eps = 1e-05


@pytest.mark.parametrize("model", list(np.arange(1, 19, 1)))
def test_3d_model(model):
    N_size = 64
    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library3D = os.path.join(path, "phantomlib", "Phantom3DLibrary.dat")

    # check parameters
    params = _check_params3d(model, path_library3D)
    phantom_3D = Model(model, N_size, path_library3D)

    assert np.sum(params) == 12
    assert params.dtype == np.int32
    assert 0.001 <= np.max(phantom_3D) <= 100
    assert phantom_3D.dtype == np.float32
    assert phantom_3D.shape == (N_size, N_size, N_size)


def test_3d_model_sm():
    phantom_size = (64, 32, 64)
    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library3D = os.path.join(path, "phantomlib", "Phantom3DLibrary.dat")

    # check parameters
    params = _check_params3d(11, path_library3D)
    phantom_3D = Model(11, phantom_size, path_library3D)

    assert np.sum(params) == 12
    assert params.dtype == np.int32
    assert 0.001 <= np.max(phantom_3D) <= 100
    assert phantom_3D.dtype == np.float32
    assert phantom_3D.shape == phantom_size


def test_3d_model_sub():
    phantom_size = (64, 64, 64)
    subset_ind = (15, 30)
    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library3D = os.path.join(path, "phantomlib", "Phantom3DLibrary.dat")

    # check parameters
    params = _check_params3d(11, path_library3D)
    phantom_3D = ModelSub(11, phantom_size, subset_ind, path_library3D)

    assert np.sum(params) == 12
    assert params.dtype == np.int32
    assert 0.001 <= np.max(phantom_3D) <= 100
    assert phantom_3D.dtype == np.float32
    assert phantom_3D.shape == (15, 64, 64)


@pytest.mark.parametrize("model", list(np.arange(100, 103, 1)))
def test_3d_model_temporal(model):
    N_size = 64
    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library3D = os.path.join(path, "phantomlib", "Phantom3DLibrary.dat")

    # check parameters
    params = _check_params3d(model, path_library3D)
    phantom_4D = ModelTemporal(model, N_size, path_library3D)

    assert len(params) == 12
    assert params.dtype == np.int32
    assert 0.001 <= np.max(phantom_4D) <= 100
    assert phantom_4D.dtype == np.float32
    assert phantom_4D.shape == (params[3], N_size, N_size, N_size)


def test_3d_temp_model_sub():
    phantom_size = (64, 64, 64)
    subset_ind = (15, 30)
    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library3D = os.path.join(path, "phantomlib", "Phantom3DLibrary.dat")

    # check parameters
    params = _check_params3d(100, path_library3D)
    phantom_3D = ModelTemporalSub(100, phantom_size, subset_ind, path_library3D)

    assert len(params) == 12
    assert params.dtype == np.int32
    assert 0.001 <= np.max(phantom_3D) <= 100
    assert phantom_3D.dtype == np.float32
    assert phantom_3D.shape == (params[3], 15, 64, 64)

import numpy as np
from numpy.testing import assert_allclose

import pytest
import os
import tomophantom
from tomophantom.TomoP3D import (
    _check_params3d,
    ModelSino,
    ModelSinoSub,
    ModelSinoTemporal,
    ModelSinoTemporalSub,
)

eps = 1e-05


@pytest.mark.parametrize("model", list(np.arange(1, 19, 1)))
def test_3d_models_sino(model):
    N_size = 64  # set the desired dimension of the phantom
    # Projection geometry related parameters:
    Horiz_det = int(np.sqrt(2) * N_size)  # detector column count (horizontal)
    Vert_det = N_size  # detector row count (vertical) (no reason for it to be > N)
    angles_num = int(0.5 * np.pi * N_size)
    # angles number
    angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")  # in degrees

    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library3D = os.path.join(path, "phantomlib", "Phantom3DLibrary.dat")

    # Generatea 3d sino
    sino3d = ModelSino(model, N_size, Horiz_det, Vert_det, angles, path_library3D)
    assert 0.0 <= np.max(sino3d) <= 1000
    assert sino3d.dtype == np.float32
    assert sino3d.shape == (Vert_det, angles_num, Horiz_det)


def test_3d_models_sino_sub():
    N_size = 64  # set the desired dimension of the phantom
    # Projection geometry related parameters:
    Horiz_det = int(np.sqrt(2) * N_size)  # detector column count (horizontal)
    Vert_det = N_size  # detector row count (vertical) (no reason for it to be > N)
    angles_num = int(0.5 * np.pi * N_size)
    # angles number
    angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")  # in degrees

    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library3D = os.path.join(path, "phantomlib", "Phantom3DLibrary.dat")
    subset_tuple = (15, 30)

    # Generate a 3d sino
    sino3d = ModelSinoSub(
        11, N_size, Horiz_det, Vert_det, subset_tuple, angles, path_library3D
    )
    assert 0.0 <= np.max(sino3d) <= 1000
    assert sino3d.dtype == np.float32
    assert sino3d.shape == (15, angles_num, Horiz_det)


@pytest.mark.parametrize("model", list(np.arange(100, 103, 1)))
def test_3d_models_sino_temporal(model):
    N_size = 64  # set the desired dimension of the phantom
    # Projection geometry related parameters:
    Horiz_det = int(np.sqrt(2) * N_size)  # detector column count (horizontal)
    Vert_det = N_size  # detector row count (vertical) (no reason for it to be > N)
    angles_num = int(0.5 * np.pi * N_size)
    # angles number
    angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")  # in degrees

    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library3D = os.path.join(path, "phantomlib", "Phantom3DLibrary.dat")

    params = _check_params3d(model, path_library3D)

    # Generatea 4d sino
    sino4d = ModelSinoTemporal(
        model, N_size, Horiz_det, Vert_det, angles, path_library3D
    )
    assert 0.0 <= np.max(sino4d) <= 1000
    assert sino4d.dtype == np.float32
    assert sino4d.shape == (params[3], Vert_det, angles_num, Horiz_det)


def test_3d_models_sino_temporal_sub():
    N_size = 64  # set the desired dimension of the phantom
    model = 100
    # Projection geometry related parameters:
    Horiz_det = int(np.sqrt(2) * N_size)  # detector column count (horizontal)
    Vert_det = N_size  # detector row count (vertical) (no reason for it to be > N)
    angles_num = int(0.5 * np.pi * N_size)
    # angles number
    angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")  # in degrees

    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library3D = os.path.join(path, "phantomlib", "Phantom3DLibrary.dat")
    params = _check_params3d(model, path_library3D)
    subset_tuple = (15, 30)

    # Generate a 3d sino
    sino4d = ModelSinoTemporalSub(
        model, N_size, Horiz_det, Vert_det, subset_tuple, angles, path_library3D
    )
    assert 0.0 <= np.max(sino4d) <= 1000
    assert sino4d.dtype == np.float32
    assert sino4d.shape == (params[3], 15, angles_num, Horiz_det)

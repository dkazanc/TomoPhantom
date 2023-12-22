# Defines common fixtures and makes them available to all tests

import os
import numpy as np
import pytest
import tomophantom


@pytest.fixture(scope="session")
def sino_model11_64():
    from tomophantom.TomoP2D import ModelSino

    N_size = 64  # set the desired dimension of the phantom
    model = 11
    angles_num = int(0.5 * np.pi * N_size)
    # angles number
    angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")
    P = int(np.sqrt(2) * N_size)  # detectors

    # one can specify an exact path to the parameters file
    path = os.path.dirname(tomophantom.__file__)
    path_library2D = os.path.join(path, "phantomlib", "Phantom2DLibrary.dat")

    # Generate sinogram
    sino_an = ModelSino(model, N_size, P, angles, path_library2D)

    return (sino_an, angles_num, P)


@pytest.fixture(scope="session")
def sino_model11_3d_64():
    from tomophantom.TomoP3D import ModelSino

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
    sino3d = ModelSino(11, N_size, Horiz_det, Vert_det, angles, path_library3D)

    return (sino3d, Vert_det, angles_num, Horiz_det)

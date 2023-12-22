import numpy as np
from numpy.testing import assert_allclose

from tomophantom.generator import foam2D, foam3D

eps = 1e-05


def test_2d_generator():
    N_size = 64  # set the desired dimension of the phantom
    tot_objects = 20  # the total number of objects to generate

    # define ranges for parameters
    x0min = -0.9
    x0max = 0.9
    y0min = -0.9
    y0max = 0.9
    c0min = 0.01
    c0max = 1.0
    ab_min = 0.01
    ab_max = 0.25

    # 2D example
    (Objfoam2D, myObjects) = foam2D(
        x0min,
        x0max,
        y0min,
        y0max,
        c0min,
        c0max,
        ab_min,
        ab_max,
        N_size,
        tot_objects,
        object_type="mix",
    )

    assert 0.0 <= np.max(Objfoam2D) <= 10
    assert Objfoam2D.dtype == np.float32
    assert Objfoam2D.shape == (N_size, N_size)


def test_3d_generator():
    N_size = 64  # define the grid size
    tot_objects = 10  # the total number of objects to generate

    # define ranges for parameters
    x0min = -0.9
    x0max = 0.9
    y0min = -0.9
    y0max = 0.9
    z0min = -0.9
    z0max = 0.9
    c0min = 0.01
    c0max = 1.0
    ab_min = 0.01
    ab_max = 0.25

    (Objfoam3D, myObjects) = foam3D(
        x0min,
        x0max,
        y0min,
        y0max,
        z0min,
        z0max,
        c0min,
        c0max,
        ab_min,
        ab_max,
        N_size,
        tot_objects,
        object_type="mix",
    )

    assert 0.0 <= np.max(Objfoam3D) <= 10
    assert Objfoam3D.dtype == np.float32
    assert Objfoam3D.shape == (N_size, N_size, N_size)

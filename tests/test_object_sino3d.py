import numpy as np
from numpy.testing import assert_allclose

from tomophantom.TomoP3D import Objects3D, ObjectSino


def test_3d_object_sino():
    N3D_size = 64  # set the desired dimension of the phantom
    Horiz_det = int(np.sqrt(2) * N3D_size)  # detector column count (horizontal)
    Vert_det = N3D_size  # detector row count (vertical) (no reason for it to be > N)
    angles_num = int(0.5 * np.pi * N3D_size)
    # angles number
    angles = np.linspace(0.0, 179.9, angles_num, dtype="float32")  # in degrees

    # specify object parameters, here we replicate model
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

    myObjects = [obj3D_1, obj3D_2]  # dictionary of objects
    ProjData3D = ObjectSino(N3D_size, Horiz_det, Vert_det, angles, myObjects)

    assert 0.0 <= np.max(ProjData3D) <= 1000
    assert ProjData3D.dtype == np.float32
    assert ProjData3D.shape == (Vert_det, angles_num, Horiz_det)

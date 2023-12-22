import numpy as np

from tomophantom.flatsgen import synth_flats


def test_3d_models_sino(sino_model11_3d_64):
    sino3d, Vert_det, angles_num, Horiz_det = sino_model11_3d_64

    I0 = 50000
    # full-beam photon flux intensity
    flatsnum = 20  # the number of the flat fields required

    [projData3D_raw, flats_combined3D, speckel_map] = synth_flats(
        sino3d,
        source_intensity=I0,
        detectors_miscallibration=0.05,
        variations_number=3,
        arguments_Bessel=(1, 25),
        specklesize=2,
        kbar=2,
        jitter_projections=0.0,
        sigmasmooth=3,
        flatsnum=flatsnum,
    )

    assert 0.0 <= np.max(projData3D_raw) <= 65535
    assert projData3D_raw.dtype == np.uint16
    assert projData3D_raw.shape == (Vert_det, angles_num, Horiz_det)
    assert flats_combined3D.dtype == np.uint16
    assert flats_combined3D.shape == (Vert_det, flatsnum, Horiz_det)

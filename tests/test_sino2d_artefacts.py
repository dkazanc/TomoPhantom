import numpy as np
from numpy.testing import assert_allclose

import pytest
from tomophantom.artefacts import artefacts_mix

eps = 1e-05


@pytest.mark.parametrize("noise_type", ["Poisson", "Gaussian"])
def test_2d_sino_noise(noise_type, sino_model11_64):
    sino_an, angles_num, P = sino_model11_64
    ideal_sino_sum = np.sum(sino_an)
    if noise_type == "Poisson":
        _noise_ = {
            "noise_type": noise_type,
            "noise_amplitude": 10000,
            "noise_seed": 1,
        }
    else:
        _noise_ = {
            "noise_type": noise_type,
            "noise_amplitude": 0.1,
            "noise_seed": 1,
        }
    noisy_sino = artefacts_mix(sino_an, **_noise_)
    noisy_sino_sum = np.sum(noisy_sino)
    assert noisy_sino_sum != ideal_sino_sum
    assert 0.0 <= np.max(noisy_sino) <= 500
    assert noisy_sino.dtype == np.float32
    assert noisy_sino.shape == (angles_num, P)


def test_2d_sino_zingers(sino_model11_64):
    sino_an, angles_num, P = sino_model11_64
    ideal_sino_sum = np.sum(sino_an)

    _zingers_ = {"zingers_percentage": 2, "zingers_modulus": 10}

    noisy_sino = artefacts_mix(sino_an, **_zingers_)
    noisy_sino_sum = np.sum(noisy_sino)
    assert noisy_sino_sum != ideal_sino_sum
    assert 0.0 <= np.max(noisy_sino) <= 500
    assert noisy_sino.dtype == np.float32
    assert noisy_sino.shape == (angles_num, P)


def test_2d_sino_stripes(sino_model11_64):
    sino_an, angles_num, P = sino_model11_64
    ideal_sino_sum = np.mean(sino_an)

    _stripes_ = {
        "stripes_percentage": 0.8,
        "stripes_maxthickness": 2.0,
        "stripes_intensity": 0.25,
        "stripes_type": "full",
        "stripes_variability": 0.002,
    }

    noisy_sino = artefacts_mix(sino_an, **_stripes_)
    noisy_sino_sum = np.sum(np.abs(noisy_sino - ideal_sino_sum))
    assert noisy_sino_sum > 0.0
    assert 0.0 <= np.max(noisy_sino) <= 500
    assert noisy_sino.dtype == np.float32
    assert noisy_sino.shape == (angles_num, P)


def test_2d_sino_shifts_pixel(sino_model11_64):
    sino_an, angles_num, P = sino_model11_64

    _sinoshifts_ = {"datashifts_maxamplitude_pixel": 10}

    [noisy_sino_misalign, shifts] = artefacts_mix(sino_an, **_sinoshifts_)
    noisy_sino_sum = np.sum(np.abs(noisy_sino_misalign - sino_an))
    assert noisy_sino_sum > 0.0
    assert 0.0 <= np.max(noisy_sino_misalign) <= 500
    assert noisy_sino_misalign.dtype == np.float32
    assert noisy_sino_misalign.shape == (angles_num, P)


def test_2d_sino_shifts_subpixel(sino_model11_64):
    sino_an, angles_num, P = sino_model11_64

    _sinoshifts_ = {"datashifts_maxamplitude_subpixel": 0.7}

    [noisy_sino_misalign, shifts] = artefacts_mix(sino_an, **_sinoshifts_)
    noisy_sino_sum = np.sum(np.abs(noisy_sino_misalign - sino_an))
    assert noisy_sino_sum > 0.0
    assert 0.0 <= np.max(noisy_sino_misalign) <= 500
    assert noisy_sino_misalign.dtype == np.float32
    assert noisy_sino_misalign.shape == (angles_num, P)


def test_2d_sino_pve(sino_model11_64):
    sino_an, angles_num, P = sino_model11_64

    _pve_ = {"pve_strength": 3}

    noisy_sino = artefacts_mix(sino_an, **_pve_)
    noisy_sino_sum = np.sum(np.abs(noisy_sino - sino_an))
    assert noisy_sino_sum > 0.0
    assert 0.0 <= np.max(noisy_sino) <= 500
    assert noisy_sino.dtype == np.float32
    assert noisy_sino.shape == (angles_num, P)


def test_2d_sino_frenel(sino_model11_64):
    sino_an, angles_num, P = sino_model11_64

    _fresnel_ = {
        "fresnel_dist_observation": 50,
        "fresnel_scale_factor": 10,
        "fresnel_wavelenght": 0.001,
    }

    noisy_sino = artefacts_mix(sino_an, **_fresnel_)
    noisy_sino_sum = np.sum(np.abs(noisy_sino - sino_an))
    assert noisy_sino_sum > 0.0
    assert 0.0 <= np.max(noisy_sino) <= 500
    assert noisy_sino.dtype == np.float32
    assert noisy_sino.shape == (angles_num, P)


def test_2d_sino_combi(sino_model11_64):
    sino_an, angles_num, P = sino_model11_64

    _noise_ = {
        "noise_type": "Poisson",
        "noise_amplitude": 10000,
        "noise_seed": 1,
    }
    _zingers_ = {"zingers_percentage": 2, "zingers_modulus": 10}

    _stripes_ = {
        "stripes_percentage": 0.8,
        "stripes_maxthickness": 2.0,
        "stripes_intensity": 0.25,
        "stripes_type": "full",
        "stripes_variability": 0.002,
    }

    _sinoshifts_ = {"datashifts_maxamplitude_subpixel": 0.7}

    _pve_ = {"pve_strength": 3}

    _fresnel_ = {
        "fresnel_dist_observation": 50,
        "fresnel_scale_factor": 10,
        "fresnel_wavelenght": 0.001,
    }

    [noisy_sino, shifts] = artefacts_mix(
        sino_an,
        **_noise_,
        **_zingers_,
        **_stripes_,
        **_sinoshifts_,
        **_pve_,
        **_fresnel_,
    )

    noisy_sino_sum = np.sum(np.abs(noisy_sino - sino_an))
    assert noisy_sino_sum > 0.0
    assert 0.0 <= np.max(noisy_sino) <= 500
    assert noisy_sino.dtype == np.float32
    assert noisy_sino.shape == (angles_num, P)

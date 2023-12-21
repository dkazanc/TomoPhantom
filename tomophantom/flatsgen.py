#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright 2023
The University of Manchester & Diamond Light Source
Licensed under the Apache License, Version 2.0.
"""

from scipy.special import spherical_yn
from scipy.ndimage import gaussian_filter
from scipy.ndimage import shift
import random
import numpy as np
from tomophantom.artefacts import noise
from tomophantom.supp.speckle_routines import simulate_speckles_with_shot_noise


def synth_flats(
    projData3D_clean: np.ndarray,
    source_intensity: int,
    detectors_miscallibration: float = 0.05,
    variations_number: int = 3,
    arguments_Bessel: tuple = (1, 25),
    specklesize: int = 2,
    kbar: int = 2,
    sigmasmooth: int = 3,
    jitter_projections: float = 0.0,
    flatsnum: int = 20,
) -> list:
    """Function to generate synthetic flat field images and raw data for projection data normalisation.
    This is more realistic modelling of stripes and various normalisation artefacts.

    Args:
        projData3D_clean (np.ndarray):  3D projection data of the following shape: (DetectorsVert, anglesDim, DetectorsHoriz).
        source_intensity (int): Source intensity which affects the amount of the Poisson noise added to data.
        detectors_miscallibration (float, optional): A constant which perturbs some detectors positions, leading to more prounonced ring artefacts. Defaults to 0.05.
        variations_number (int, optional): A type of function to control stripe type (1 - linear, 2, - sinusoidal, 3 - exponential). Defaults to 3.
        arguments_Bessel (tuple, optional): A tuple of 2 arguments for 2 Bessel functions to control background variations in flats. Defaults to (1,25).
        specklesize (int, optional):  Speckle size in pixel units for flats background simulation. Defaults to 2.
        kbar (int, optional): Mean photon density (photons per pixel) for background simulation. Defaults to 2.
        sigmasmooth (int, optional): Gaussian smoothing parameter to blur the speckled backround (1,3,5,7,...). Defaults to 3.
        jitter_projections (float, optional): An additional random jitter to the projections in pixels. Defaults to 0.0.
        flatsnum (int, optional): A number of flats to generate. Defaults to 20.

    Returns:
        list: List that contains [np.uint16(projData3D_raw), np.uint16(simulated_flats), blurred_speckles_map]
    """
    [DetectorsDimV, projectionsNo, DetectorsDimH] = np.shape(projData3D_clean)

    # output datasets
    flats_combined3D = np.zeros(
        (DetectorsDimV, flatsnum, DetectorsDimH), dtype="uint16"
    )
    projData3D_raw = np.zeros(np.shape(projData3D_clean), dtype="float32")

    # normalise the data
    projData3D_clean /= np.max(projData3D_clean)

    # using spherical Bessel functions to emulate the background (scintillator) variations
    func = spherical_yn(
        1,
        np.linspace(
            arguments_Bessel[0], arguments_Bessel[1], DetectorsDimV, dtype="float32"
        ),
    )
    func += abs(np.min(func))

    flatfield = np.zeros((DetectorsDimV, DetectorsDimH))
    for i in range(0, DetectorsDimH):
        flatfield[:, i] = func
    for i in range(0, DetectorsDimH):
        flatfield[:, i] += np.flipud(func)

    if specklesize != 0.0:
        # using speckle generator routines to create a photon count texture in the background
        speckle_background = simulate_speckles_with_shot_noise(
            [DetectorsDimV, DetectorsDimH], 1, specklesize, kbar
        )
    else:
        speckle_background = np.ones((DetectorsDimV, DetectorsDimH))

    # model miscallibrated detectors (a possible path to generate ring artifacts)
    blurred_speckles_map = np.zeros((DetectorsDimV, DetectorsDimH, variations_number))
    for i in range(0, variations_number):
        speckles = simulate_speckles_with_shot_noise(
            [DetectorsDimV, DetectorsDimH], 1, 10, 0.03
        )
        # blur the speckled background
        blurred_speckles = gaussian_filter(speckles.copy(), sigma=sigmasmooth)
        # threshold the result
        blurred_speckles[blurred_speckles < 0.6 * np.max(blurred_speckles)] = 0
        blurred_speckles_map[:, :, i] = blurred_speckles
    blurred_speckles_map /= np.max(blurred_speckles_map)

    sinusoidal_response = (
        np.sin(np.linspace(0, 1.5 * np.pi, projectionsNo))
        + np.random.random(projectionsNo) * 0.1
    )
    sinusoidal_response /= np.max(sinusoidal_response)
    exponential_response = (
        np.exp(np.linspace(0, np.pi, projectionsNo))
        + np.random.random(projectionsNo) * 0.1
    )
    exponential_response /= np.max(exponential_response)

    # prepeare flat fields
    for i in range(0, flatsnum):
        # add speckled background to the initial image with the Bessel background
        flatfield_combined = flatfield.copy() + 0.5 * (
            speckle_background / np.max(speckle_background)
        )
        flatfield_combined /= np.max(flatfield_combined)

        # adding Poisson noise to flat fields
        flatfield_poisson = noise(
            flatfield_combined * source_intensity, source_intensity, noisetype="Poisson"
        )
        flatfield_poisson /= np.max(flatfield_poisson)

        flats_combined3D[:, i, :] = np.uint16(flatfield_poisson * 65535)

    # convert synthetic projections to raw-data like projection ready for normalisation
    for i in range(0, projectionsNo):
        proj_exp = (
            np.exp(-projData3D_clean[:, i, :]) * source_intensity * flatfield_poisson
        )  # raw projection
        for j in range(0, variations_number):
            if j == 0:
                # adding a consistent offset for certain detectors
                proj_exp += (
                    blurred_speckles_map[:, :, j]
                    * detectors_miscallibration
                    * source_intensity
                )
            if j == 1:
                # adding a sinusoidal-like response offset for certain detectors
                proj_exp += (
                    sinusoidal_response[i]
                    * blurred_speckles_map[:, :, j]
                    * detectors_miscallibration
                    * source_intensity
                )
            if j == 2:
                # adding an exponential response offset for certain detectors
                proj_exp += (
                    exponential_response[i]
                    * blurred_speckles_map[:, :, j]
                    * detectors_miscallibration
                    * source_intensity
                )

        projection_poisson = noise(proj_exp, source_intensity, noisetype="Poisson")

        # apply jitter to projections
        if jitter_projections != 0.0:
            horiz_shift = random.uniform(
                -jitter_projections, jitter_projections
            )  # generate random directional shift
            vert_shift = random.uniform(
                -jitter_projections, jitter_projections
            )  # generate random directional shift
            projection_poisson = shift(
                projection_poisson.copy(), [vert_shift, horiz_shift], mode="reflect"
            )
        projData3D_raw[:, i, :] = projection_poisson

    projData3D_raw /= np.max(projData3D_raw)
    return [
        np.uint16(projData3D_raw * 65535),
        np.uint16(flats_combined3D),
        blurred_speckles_map,
    ]

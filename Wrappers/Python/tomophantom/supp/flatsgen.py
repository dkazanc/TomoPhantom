#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A function to generate synthetic flat field images for 3D projection data normalisation

@author: Daniil Kazantsev
"""

from scipy.special import spherical_yn
from scipy.special import y1
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.interpolation import shift
import random
import numpy as np
from tomophantom.supp.artifacts import noise
from tomophantom.supp.speckle_routines import simulate_speckles_with_shot_noise


def synth_flats(projData3D_clean, source_intensity, source_variation, arguments_Bessel, specklesize, kbar, sigmasmooth, jitter, flatsnum):
    """
    the required format of the input (clean) data is [detectorsX, Projections, detectorsY]
    Parameters: 
    source_intensity - source intensity which affects the amount of Poisson noise added to data
    source_variation - constant which perturbs the source intensity leading to ring artifacts etc.
    arguments_Bessel - tuple of 4 Arguments for 2 Bessel functions to control background variations
    specklesize - speckle size in pixel units for background simulation
    kbar  - mean photon density (photons per pixel) for background simulation
    jitter - a random jitter to the speckled background given in pixels
    sigmasmooth - Gaussian smoothing parameter to blur the speckled backround (1,3,5,7...)
    flatsnum - a number of flats to generate
    """
    [DetectorsDimV, projectionsNo, DetectorsDimH] = np.shape(projData3D_clean)
    flatfield = np.zeros((DetectorsDimV,DetectorsDimH))
    blurred_speckles = np.zeros((DetectorsDimV,DetectorsDimH))
    blurred_speckles_res = np.zeros((DetectorsDimV,DetectorsDimH))
    flats_combined3D = np.zeros((DetectorsDimV,flatsnum, DetectorsDimH))
    projData3D_noisy = np.zeros(np.shape(projData3D_clean),dtype='float32')
    source_intensity_var = source_intensity*np.ones((DetectorsDimV,DetectorsDimH))
    source_intensity_variable = source_intensity_var
    maxProj_scalar = np.max(projData3D_clean)
    
    # using spherical Bessel functions to emulate the background (scintillator) variations
    func = spherical_yn(1, np.linspace(arguments_Bessel[0], arguments_Bessel[1], DetectorsDimV,dtype='float32'))
    func = func + abs(np.min(func))
    func2 = y1(np.linspace(arguments_Bessel[2],arguments_Bessel[3],DetectorsDimH,dtype='float32'))
    func2 = func2 + abs(np.min(func2))
    
    for i in range(0,DetectorsDimV):
        flatfield[i,:] = func2
    
    for i in range(0,DetectorsDimH):
        flatfield[:,i] += func
    
    if (specklesize != 0.0):
        # using speckle generator routines to create a texture in the background
        modes = 1
        speckle_background = simulate_speckles_with_shot_noise([DetectorsDimV, DetectorsDimH], modes, specklesize, kbar)
        #blur the speckled background and add to the initial image with the Bessel background
        blurred_speckles = gaussian_filter(speckle_background.copy(), sigma=sigmasmooth)

    for i in range(0,flatsnum):
        # adding noise and normalise
        if (jitter != 0.0):
            horiz_shift = random.uniform(-jitter,jitter)  #generate random directional shift
            vert_shift = random.uniform(-jitter,jitter)  #generate random directional shift
            blurred_speckles_res = shift(blurred_speckles.copy(),[vert_shift,horiz_shift])
        else:
            blurred_speckles_res = blurred_speckles
        flat_combined = flatfield + blurred_speckles_res
        flat_combined /= np.max(flat_combined)
        # make source intensity variable if required
        if (source_variation is not None or source_variation != 0.0):
            source_intensity_variable = noise(source_intensity_var, source_variation*source_intensity, noisetype='Gaussian', seed = None, prelog = None)

        #adding Poisson noise to the flat fields
        flat_noisy = np.random.poisson(np.multiply(source_intensity_variable,flat_combined))
        flats_combined3D[:,i,:] = flat_noisy

    for i in range(0,projectionsNo):
        # make source intensity variable if required
        if (source_variation is not None or source_variation != 0.0):
            source_intensity_variable = noise(source_intensity_var, source_variation*source_intensity, noisetype='Gaussian', seed = None, prelog = None)
        # adding noise and normalise
        if (jitter != 0.0):
            horiz_shift = random.uniform(-jitter,jitter)  #generate random directional shift
            vert_shift = random.uniform(-jitter,jitter)  #generate random directional shift
            blurred_speckles_res = shift(blurred_speckles.copy(),[vert_shift,horiz_shift])
        else:
            blurred_speckles_res = blurred_speckles
        flat_combined = flatfield + blurred_speckles_res
        flat_combined /= np.max(flat_combined)
        projData3D_noisy[:,i,:] = np.random.poisson(np.random.poisson(np.multiply(source_intensity_variable,flat_combined))* np.exp(-projData3D_clean[:,i,:]/maxProj_scalar))

    return [projData3D_noisy, flats_combined3D]
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Daniil Kazantsev
"""

from scipy.special import spherical_yn
from scipy.ndimage import gaussian_filter
from scipy.ndimage import shift
import random
import numpy as np
from tomophantom.supp.artifacts import noise
from tomophantom.supp.speckle_routines import simulate_speckles_with_shot_noise


def synth_flats(projData3D_clean, source_intensity, detectors_miscallibration = 0.05, variations_number = 3, arguments_Bessel=(1,25), specklesize = 2, kbar = 2, sigmasmooth = 2, jitter_projections = 0.0, flatsnum=20):
    """
    A function to generate synthetic flat field images and raw data for projection data normalisation.
    This is a way to more realistic modelling of stripes leading to ring artifacts
    the format of the input (clean) data is [detectorsX, Projections, detectorsY]
    
    Parameters:
    source_intensity - source intensity which affects the amount of the Poisson noise added to data
    variations_number - the number of functions to control stripe type (1 - linear, 2, - sinusoidal, 3 - exponential)
    detectors_miscallibration - a constant which perturbs the some detectors positions, leading to ring artifacts etc.
    arguments_Bessel - tuple of 2 Arguments for 2 Bessel functions to control background variations
    specklesize - speckle size in pixel units for background simulation
    kbar  - mean photon density (photons per pixel) for background simulation
    jitter_projections - a random jitter to the projections in pixels
    sigmasmooth - Gaussian smoothing parameter to blur the speckled backround (1,3,5,7...)
    flatsnum - a number of flats to generate
    """
    [DetectorsDimV, projectionsNo, DetectorsDimH] = np.shape(projData3D_clean)
    
    # output datasets
    flats_combined3D = np.zeros((DetectorsDimV,flatsnum, DetectorsDimH),dtype='uint16')
    projData3D_raw = np.zeros(np.shape(projData3D_clean),dtype='float32')
    
    # normalise the data
    projData3D_clean /= np.max(projData3D_clean)
    
    # using spherical Bessel functions to emulate the background (scintillator) variations
    func = spherical_yn(1, np.linspace(arguments_Bessel[0], arguments_Bessel[1], DetectorsDimV,dtype='float32'))
    func += abs(np.min(func))
    
    flatfield = np.zeros((DetectorsDimV,DetectorsDimH))
    for i in range(0,DetectorsDimH):
        flatfield[:,i] = func
    for i in range(0,DetectorsDimH):
        flatfield[:,i] += np.flipud(func)
    
    if (specklesize != 0.0):
        # using speckle generator routines to create a photon count texture in the background
        speckle_background = simulate_speckles_with_shot_noise([DetectorsDimV, DetectorsDimH], 1, specklesize, kbar)
    else:
        speckle_background = np.ones((DetectorsDimV,DetectorsDimH))

    # model miscallibrated detectors (a possible path to generate ring artifacts)
    blurred_speckles_map = np.zeros((DetectorsDimV, DetectorsDimH, variations_number))
    for i in range(0,variations_number):
        speckles = simulate_speckles_with_shot_noise([DetectorsDimV, DetectorsDimH], 1, 10, 0.03)
        #blur the speckled background
        blurred_speckles = gaussian_filter(speckles.copy(), sigma=sigmasmooth)
        # threshold the result
        blurred_speckles[blurred_speckles < 0.6*np.max(blurred_speckles)] = 0
        blurred_speckles_map[:,:,i] = blurred_speckles
    blurred_speckles_map /= np.max(blurred_speckles_map)
    
    sinusoidal_response = np.sin(np.linspace(0,1.5*np.pi,projectionsNo)) + np.random.random(projectionsNo) * 0.1
    sinusoidal_response /= np.max(sinusoidal_response)
    exponential_response = np.exp(np.linspace(0,np.pi,projectionsNo)) + np.random.random(projectionsNo) * 0.1
    exponential_response /= np.max(exponential_response)

    # prepeare flat fields
    for i in range(0,flatsnum):
        # add speckled background to the initial image with the Bessel background
        flatfield_combined = flatfield.copy() + 0.5*(speckle_background/np.max(speckle_background))
        flatfield_combined /= np.max(flatfield_combined)
        
        #adding Poisson noise to flat fields
        flatfield_poisson = noise(flatfield_combined*source_intensity, source_intensity, noisetype='Poisson')
        flatfield_poisson /= np.max(flatfield_poisson)

        flats_combined3D[:,i,:] = np.uint16(flatfield_poisson*65535)

    # convert synthetic projections to raw-data like projection ready for normalisation
    for i in range(0,projectionsNo):
        proj_exp = np.exp(-projData3D_clean[:,i,:])*source_intensity*flatfield_poisson # raw projection
        for j in range(0,variations_number):
            if j == 0:
                # adding a consistent offset for certain detectors
                proj_exp += blurred_speckles_map[:,:,j]*detectors_miscallibration*source_intensity
            if j == 1:
                # adding a sinusoidal-like response offset for certain detectors
                proj_exp += sinusoidal_response[i]*blurred_speckles_map[:,:,j]*detectors_miscallibration*source_intensity
            if j == 2:
                # adding an exponential response offset for certain detectors
                proj_exp += exponential_response[i]*blurred_speckles_map[:,:,j]*detectors_miscallibration*source_intensity
                
        projection_poisson = noise(proj_exp, source_intensity, noisetype='Poisson')

        # apply jitter to projections
        if jitter_projections != 0.0:
            horiz_shift = random.uniform(-jitter_projections,jitter_projections)  #generate random directional shift
            vert_shift = random.uniform(-jitter_projections,jitter_projections)  #generate random directional shift
            projection_poisson = shift(projection_poisson.copy(),[vert_shift,horiz_shift], mode='reflect')
        projData3D_raw[:,i,:] = projection_poisson

    projData3D_raw /= np.max(projData3D_raw)
    return [np.uint16(projData3D_raw*65535), np.uint16(flats_combined3D), blurred_speckles_map]
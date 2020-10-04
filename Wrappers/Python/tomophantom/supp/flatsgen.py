#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A function to generate synthetic flat field images for 3D projection data normalisation

@author: Daniil Kazantsev
"""

from scipy.special import spherical_yn
from scipy.special import y1
from scipy.ndimage.filters import gaussian_filter
import random
import numpy as np
from tomophantom.supp.artifacts import noise

def synth_flats(projData3D_clean, source_intensity, source_variation, arguments_Bessel, strip_height, strip_thickness, sigmasmooth, flatsnum):
    """
    the required format of the input (clean) data is [detectorsX, Projections, detectorsY]
    Parameters: 
    source_intensity - a source intensity which affects the amount of Poisson noise added to data
    source_variation - a constant which perturbs the source intensity leading to ring artifacts etc.
    arguments_Bessel - a tuple of 4 Arguments for 2 Bessel functions to control background variations
    strip_height - a value between (0,1] which controls the height of the background stripes
    strip_thickness - a value in pixels which controls the width of stripes
    sigmasmooth - a smoothing parameter for a stripe image with noise (1,3,5,7...)
    flatsnum - a number of flats to generate
    """
    [DetectorsDimV, projectionsNo, DetectorsDimH] = np.shape(projData3D_clean)
    flatfield = np.zeros((DetectorsDimV,DetectorsDimH))
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
    
    # Adding stripes of varying vertical & horizontal sizes
    maxheight_par = round(strip_height*DetectorsDimV)
    max_intensity = np.max(flatfield)
    flatstripes = np.zeros(np.shape(flatfield))
    
    for x in range(0,DetectorsDimH):
        randind = random.randint(0,DetectorsDimV) # generate random index (vertically)
        randthickness = random.randint(0,strip_thickness) #generate random thickness
        randheight = random.randint(0,maxheight_par) #generate random height
        randintens = random.uniform(0.1, 2.0) # generate random multiplier
        intensity = max_intensity*randintens
        
        for x1 in range(-randheight,randheight):
            if (((randind+x1) >= 0) and ((randind+x1) < DetectorsDimV)):
                for x2 in range(-randthickness,randthickness):
                    if (((x+x2) >= 0) and ((x+x2) < DetectorsDimH)):
                        flatstripes[randind+x1,x+x2] += intensity
    
    for i in range(0,flatsnum):
        # adding noise and normalise
        flatstripes2 = flatstripes.copy()
        #blur the result and add to the initial image with the Bessel background
        blurred_flatstripes = gaussian_filter(flatstripes2, sigma=sigmasmooth)
        flat_combined = flatfield + blurred_flatstripes
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
        projData3D_noisy[:,i,:] = np.random.poisson(np.random.poisson(np.multiply(source_intensity_variable,flat_combined))* np.exp(-projData3D_clean[:,i,:]/maxProj_scalar))

    return [projData3D_noisy, flats_combined3D]
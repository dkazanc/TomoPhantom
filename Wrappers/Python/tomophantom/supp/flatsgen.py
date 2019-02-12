#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A function to generate fake flat field images for 3D projection data normalisation

@author: Daniil Kazantsev
"""

from scipy.special import spherical_yn
from scipy.special import y1
from scipy.ndimage.filters import gaussian_filter
import random
import numpy as np
# import matplotlib.pyplot as plt

def flats(DetectorsDimV, DetectorsDimH, maxheight, maxthickness, sigma_noise, sigmasmooth, flatsnum):
    """
    maxheight - a value between (0,1] which controls the height of stripes
    maxthickness - a value in pixels which controls the width of stripes
    sigma_noise - a noise level (Gaussian) which is added to the sripe image
    sigmasmooth - a smoothing parameter for a stripe image with noise (1,3,5,7...)
    flatsnum - a number of flats to generate
    """
    
    flatfield = np.zeros((DetectorsDimV,DetectorsDimH))
    flat_combined3D = np.zeros((flatsnum,DetectorsDimV,DetectorsDimH))
    
    # using spherical Bessel functions
    func = spherical_yn(1, np.linspace(3,15,DetectorsDimV,dtype='float32'))
    func = func + abs(np.min(func))
    func2 = y1(np.linspace(15,5,DetectorsDimH,dtype='float32'))
    func2 = func2 + abs(np.min(func2))
    
    for i in range(0,DetectorsDimV):
        flatfield[i,:] = func2
    
    for i in range(0,DetectorsDimH):
        flatfield[:,i] += func
    
    # Adding stripes of varying vertical & horizontal sizes
    maxheight_par = round(maxheight*DetectorsDimV)
    max_intensity = np.max(flatfield)
    flatstripes = np.zeros(np.shape(flatfield))
    
    for x in range(0,DetectorsDimH):
        randind = random.randint(0,DetectorsDimV) # generate random index (vertically)
        randthickness = random.randint(0,maxthickness) #generate random thickness
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
        flatstripes2 += np.random.normal(loc = 0.00 ,scale = sigma_noise, size = np.shape(flatstripes2))
        flatstripes2 /= np.max(flatstripes2)
        flatstripes2 += abs(np.min(flatstripes2))
        
        #blur the result and add to the initial image with bessel background
        blurred_flatstripes = gaussian_filter(flatstripes2, sigma=sigmasmooth)
        flat_combined = flatfield + blurred_flatstripes
        flat_combined /= np.max(flat_combined)
        flat_combined3D[i,:,:] = flat_combined
    return flat_combined3D
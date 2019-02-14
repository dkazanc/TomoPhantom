#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A function to normalise projection data using the generated flat fields, Poisson
noise is applied to flat fields and data

The errors in the flat fields will be propagated into prokjection data after
normalisation. This can have a serious negative effect on the reconstructed images 
in terms of the quality of reconstruction (artifacts presence). 

This modelling, however, is more realistic and recommended if algorithms are 
tested on robustness towards artifacts and errors in data

@author: Daniil Kazantsev
"""
import random
import numpy as np
#import matplotlib.pyplot as plt
from tomophantom.supp.artifacts import ArtifactsClass

def normaliser_sim(projData3D, flatSIM, sigma_flats = 0.05, flux_intensity = 30000):
    """
    projData3D - 3D projection data (noiseless) [DetectorsDimV, Proj_angles, DetectorsDimH]
    maxthickness - a value in pixels which controls the width of stripes
    sigma_flats - a noise level (Gaussian) in flats, do not set too high to avoid outliers
    flux_intensity  - controls the level of Posson noise applied to projection data
    """
    [DetectorsDimV, Proj_angles, DetectorsDimH] = np.shape(projData3D)
    [flatsnum, DetectorsDimV_f, DetectorsDimH_f] = np.shape(flatSIM)
    if (DetectorsDimV != DetectorsDimV_f):
        raise("The size of the vertical detector for data and the flat field is different ")
    if (DetectorsDimH != DetectorsDimH_f):
        raise("The size of the horizontal detector for data and the flat field is different ")
    # add noise to the stack of flat images
    artifacts_add = ArtifactsClass(flatSIM)
    flat_all_noise = artifacts_add.noise(sigma=sigma_flats,noisetype='Gaussian')
    flat_average_noise = np.average(flat_all_noise,0) # calculate average of all flats

    nonzeroInd = np.where(flat_average_noise != 0) # nonzero data
    zeroInd = np.where(flat_average_noise == 0) # zero data
    projData3D_norm = np.zeros(np.shape(projData3D),dtype='float32')
    
    for x in range(0,Proj_angles):
        proj2D = projData3D[:,x,:]
        norm_proj = np.zeros(np.shape(proj2D))
        randflatind = random.randint(0,flatsnum-1) # generate random flat number
        proj_flat = flatSIM[randflatind,:,:]*proj2D
        artifacts_add = ArtifactsClass(proj_flat) # adding Poisson noise
        proj_flat_noisy = artifacts_add.noise(sigma=flux_intensity,noisetype='Poisson')
        norm_proj[nonzeroInd]  = proj_flat_noisy[nonzeroInd]/flat_average_noise[nonzeroInd]
        norm_proj[zeroInd] = 1e-13
        projData3D_norm[:,x,:] = np.float32(norm_proj)

    return projData3D_norm
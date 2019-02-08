#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A function to normalise projection data using the generated flat fields. 
Poisson noise is applied to flat fields and data
"""
import random
import numpy as np
import matplotlib.pyplot as plt
from tomophantom.supp.artifacts import ArtifactsClass
from flatsgen import flats

# def normaliser2(DetectorsDimV, DetectorsDimH, maxheight, maxthickness, sigma_noise, sigmasmooth, flatsnum):
"""
maxheight - a value between (0,1] which controls the height of stripes

"""
flux_intensity = 30000
flatsnum = 20
sliceSel = 128

[DetectorsDimV, DetectorsDimH] = np.shape(proj2D)
flatsynth = flats(DetectorsDimV, DetectorsDimH, maxheight = 0.1, maxthickness = 3, sigma_noise = 0.3, sigmasmooth = 3, flatsnum=flatsnum)
flat_average = np.average(flatsynth,0)

artifacts_add = ArtifactsClass(flat_average)
flat_average_noise = artifacts_add.noise(sigma=0.05,noisetype='Gaussian')
# flat_average_noise /= np.max(flat_average_noise)
#flat_average_noise = flat_average.copy()

nonzeroInd = np.where(flat_average_noise != 0) # nonzero data
zeroInd = np.where(flat_average_noise == 0) # zero data
projData3D_norm = np.zeros(np.shape(projData3D_analyt),dtype='float32')

for x in range(0,proj_num):
    proj2D = projData3D_analyt[:,x,:]
    norm_proj = np.zeros(np.shape(proj2D))
    proj_noisy_flat = np.zeros(np.shape(proj2D))
    randflatind = random.randint(0,flatsnum-1) # generate random flat number
    selected_flat = flatsynth[randflatind,:,:]
    artifacts_add = ArtifactsClass(proj2D)
    proj2D_noisy = artifacts_add.noise(sigma=flux_intensity,noisetype='Poisson')
    proj_noisy_flat[nonzeroInd] = selected_flat[nonzeroInd]*proj2D_noisy[nonzeroInd]
    norm_proj[nonzeroInd]  = proj_noisy_flat[nonzeroInd]/flat_average_noise[nonzeroInd]
    norm_proj[zeroInd] = 1e-13
    projData3D_norm[:,x,:] = np.float32(norm_proj)

#return flat_combined3D
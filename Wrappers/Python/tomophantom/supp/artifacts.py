#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  9 10:47:01 2018
Note that the TomoPhantom package is released under Apache License, Version 2.0

Artifacts simulation class for sinograms (2D) or for 3D-projection data
 
What can be simulated: 
-- noise (Poisson or Gaussian)
-- zingers (streaks in reconstructions)
-- stripes (rings in reconstructions)
-- shifts - misalignment (blur in reconstruction)

@author: Daniil Kazantsev
"""
import numpy as np
import random

class ArtifactsClass:
    def __init__(self, sinogram):
        self.sinogram = np.copy(sinogram)
        if (self.sinogram.ndim == 2):
            (self.anglesDim, self.DetectorsDimH) = np.shape(sinogram)
        else:
            (self.DetectorsDimV, self.anglesDim, self.DetectorsDimH) = np.shape(sinogram)
    def noise(self, sigma, noisetype):
        """ adding random noise to data """
        sino_noisy = self.sinogram
        if noisetype == 'Gaussian':
            # add normal Gaussian noise
            sino_noisy += np.random.normal(loc = 0.0, scale = sigma, size = np.shape(sino_noisy))
            sino_noisy[sino_noisy<0] = 0
        elif noisetype == 'Poisson':
            # add Poisson noise
            maxSino = np.max(self.sinogram)
            if maxSino > 0:
                sino_noisy = (self.sinogram)/maxSino
                dataExp = sigma*np.exp(-sino_noisy)  # noiseless raw data
                sino_noisy = np.random.poisson(dataExp) #adding Poisson noise
                div_res = np.float32(sino_noisy)/np.max(sino_noisy)
                sino_noisy = -np.log(div_res)*maxSino #log corrected data -> sinogram
                sino_noisy[sino_noisy<0] = 0
        else:
            print ("Select 'Gaussian' or 'Poisson' for noise type")
        return sino_noisy
    def zingers(self, percentage, modulus):
        """ adding zingers (zero pixels or small 4 pixels clusters) to data 
        or 6 voxels to 3D projection data
        - percentage - the amount of zingers to be added
        - modulus controls the amount of 4/6 pixel clusters to be added
        """
        if 0 < percentage <= 100:
            pass
        else:
            raise ("percentage must be larger than zero but smaller than 100")
        if (modulus > 0):
            pass
        else:
            raise ("Modulus integer must be positive")
        sino_zingers = self.sinogram
        length_sino = np.size(sino_zingers)
        num_values = int((length_sino)*(np.float32(percentage)/100.0))
        sino_zingers_fl = sino_zingers.flatten()
        for x in range(num_values):
            randind = random.randint(0,length_sino) # generate random index 
            sino_zingers_fl[randind] =0
            if ((x % int(modulus)) == 0):
                if (self.sinogram.ndim == 2):
                    if ((randind > self.DetectorsDimH) & (randind < length_sino-self.DetectorsDimH)):
                        sino_zingers_fl[randind+1] =0
                        sino_zingers_fl[randind-1] =0
                        sino_zingers_fl[randind+self.DetectorsDimH] =0
                        sino_zingers_fl[randind-self.DetectorsDimH] =0
                else:
                    if ((randind > self.DetectorsDimH*self.DetectorsDimV) & (randind < length_sino-self.DetectorsDimH*self.DetectorsDimV)):
                        sino_zingers_fl[randind+1] =0
                        sino_zingers_fl[randind-1] =0
                        sino_zingers_fl[randind+self.DetectorsDimH] =0
                        sino_zingers_fl[randind-self.DetectorsDimH] =0
                        sino_zingers_fl[randind+self.DetectorsDimH*self.DetectorsDimV] =0
                        sino_zingers_fl[randind-self.DetectorsDimH*self.DetectorsDimV] =0
        sino_zingers[:] = sino_zingers_fl.reshape(sino_zingers.shape)
        return sino_zingers
    def stripes(self, percentage, maxthickness):
        """
        A function to add stripes (constant offsets) to sinogram which results in rings in the 
        reconstructed image
        - percentage defines the density of stripes
        - maxthickness defines maximal thickness of a stripe
        """
        if 0 < percentage <= 100:
            pass
        else:
            raise ("percentage must be larger than zero but smaller than 100")
        if 0 <= maxthickness <= 10:
            pass
        else:
            raise ("maximum thickness must be in [0,10] range")
        sino_stripes = self.sinogram
        max_intensity = np.max(sino_stripes)
        range_detect = int((np.float32(self.DetectorsDimH))*(np.float32(percentage)/100.0))
        if (self.sinogram.ndim == 2):
            for x in range(range_detect):
                randind = random.randint(0,self.DetectorsDimH) # generate random index
                randthickness = random.randint(0,maxthickness) #generate random thickness
                randintens = random.uniform(-1, 0.5) # generate random multiplier
                intensity = max_intensity*randintens
                if ((randind > 0+randthickness) & (randind < self.DetectorsDimH-randthickness)):
                    for x1 in range(-randthickness,randthickness+1):
                        sino_stripes[:,randind+x1] += intensity
        else:
            for j in range(self.DetectorsDimV):
                for x in range(range_detect):
                    randind = random.randint(0,self.DetectorsDimH) # generate random index
                    randthickness = random.randint(0,maxthickness) #generate random thickness
                    randintens = random.uniform(-1, 0.5) # generate random multiplier
                    intensity = max_intensity*randintens
                    if ((randind > 0+randthickness) & (randind < self.DetectorsDimH-randthickness)):
                        for x1 in range(-randthickness,randthickness+1):
                            sino_stripes[j,:,randind+x1] += intensity
        return sino_stripes
    def shifts(self, maxamplitude):
        """
        A function to add random shifts to sinogram rows (an offset for each angular position)
        - maxamplitude (in pixels) defines the maximal amplitude of each angular deviation
        """
        sino_shifts = np.zeros(np.shape(self.sinogram),dtype='float32')
        non = lambda s: s if s<0 else None
        mom = lambda s: max(0,s)
        for x in range(self.anglesDim):
            rand_shift = random.randint(-maxamplitude,maxamplitude)  #generate random shift
            if (self.sinogram.ndim == 2):
                projection = self.sinogram[x,:] # extract 1D projection
                projection_shift = np.zeros(np.shape(projection),dtype='float32')
                projection_shift[mom(rand_shift):non(rand_shift)] = projection[mom(-rand_shift):non(-rand_shift)]
                sino_shifts[x,:] = projection_shift
            else:
                rand_shift2 = random.randint(-maxamplitude,maxamplitude)  #generate random shift
                projection2D = self.sinogram[:,x,:] # extract 2D projection
                projection2D_shift = np.zeros(np.shape(projection2D),dtype='float32')
                projection2D_shift[mom(rand_shift):non(rand_shift), mom(rand_shift2):non(rand_shift2)] = projection2D[mom(-rand_shift):non(-rand_shift),mom(-rand_shift2):non(-rand_shift2)]
                sino_shifts[:,x,:] = projection2D_shift
        return sino_shifts

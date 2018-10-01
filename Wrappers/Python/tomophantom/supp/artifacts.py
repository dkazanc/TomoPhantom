#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  9 10:47:01 2018
Note that the TomoPhantom package is released under Apache License, Version 2.0

@author: Daniil Kazantsev

--- A class which simulates artifacts applied to sinogram (2D) data -----
 currently availble:
-- noise (Poisson or Gaussian)
-- zingers
-- stripes (rings)
"""
import numpy as np
import random

class ArtifactsClass:
    def __init__(self, sinogram):
        self.sinogram = np.copy(sinogram)
        (self.anglesDim, self.DetectorsDim) = np.shape(sinogram)
    def noise(self, sigma, noisetype):
        """ adding random noise to data """
        #sino_noisy = np.zeros(np.shape(self.sinogram),dtype='float32')
        sino_noisy = self.sinogram
        if noisetype == 'Gaussian':
            # add normal Gaussian noise
            sino_noisy += np.random.normal(loc = 0 ,scale = sigma * sino_noisy, size = np.shape(sino_noisy))
        elif noisetype == 'Poisson':
            # add Poisson noise
            maxSino = np.max(sino_noisy)
            if maxSino > 0:
                sino_noisy = sino_noisy/maxSino
                dataExp = sigma*np.exp(-sino_noisy)  # noiseless raw data
                sino_noisy = np.random.poisson(dataExp) #adding Poisson noise
                sino_noisy = -np.log(sino_noisy/np.max(sino_noisy))*maxSino #log corrected data -> sinogram
                sino_noisy[sino_noisy<0] = 0
        else:
            print ("Select 'Gaussian' or 'Poisson' for noise type")
        return sino_noisy
    def zingers(self, percentage, modulus):
        """ adding zingers (zero pixels or small 4 pixels clusters) to data 
        - percentage - the amount of zingers to be added
        - modulus controls the amount of 4 pixel clusters to be added
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
        num_values = int((length_sino)*(percentage/100))
        sino_zingers_fl = sino_zingers.flatten()
        for x in range(num_values):
            randind = random.randint(0,length_sino) # generate random index 
            sino_zingers_fl[randind] =0
            if ((x % int(modulus)) == 0):
                if ((randind > self.DetectorsDim) & (randind < length_sino-self.DetectorsDim)):
                    sino_zingers_fl[randind+1] =0
                    sino_zingers_fl[randind-1] =0
                    sino_zingers_fl[randind+self.DetectorsDim] =0
                    sino_zingers_fl[randind-self.DetectorsDim] =0
        sino_zingers[:] = sino_zingers_fl.reshape(sino_zingers.shape)
        return sino_zingers
    def stripes(self, percentage, maxthickness):
        """
        function to add stripes (constant offsets) to sinogram which results in rings in the 
        reconstructed image
        - percentage defines the density of stripes
        - maxthickness defines maximal thickness of a stripe
        """
        if 0 < percentage <= 100:
            pass
        else:
            raise ("percentage must be larger than zero but smaller than 100")
        if 0 <= maxthickness <= 5:
            pass
        else:
            raise ("maximum thickness must be in [0,5] range")
        sino_stripes = self.sinogram
        max_intensity = np.max(sino_stripes)
        range_detect = int((self.DetectorsDim)*(percentage/100))
        for x in range(range_detect):
            randind = random.randint(0,self.DetectorsDim) # generate random index
            randthickness = random.randint(0,maxthickness) #generate random thickness
            randintens = random.uniform(-1, 0.5) # generate random multiplier
            intensity = max_intensity*randintens
            if ((randind > 0+randthickness) & (randind < self.DetectorsDim-randthickness)):
                for x1 in range(-randthickness,randthickness+1):
                    sino_stripes[:,randind+x1] += intensity
        return sino_stripes


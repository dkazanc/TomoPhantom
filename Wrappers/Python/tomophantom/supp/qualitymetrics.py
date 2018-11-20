#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
A class for some standard image quality metrics

@author: Daniil Kazantsev
"""
import numpy as np

class QualityTools:
    def __init__(self, im1, im2):
        if im1.size != im2.size:
            print ('Error: Sizes of images/volumes are different')
            raise SystemExit
        self.im1 = im1 # image or volume - 1
        self.im2 = im2 # image or volume - 2
    def nrmse(self):
        """ Normalised Root Mean Square Error """
        rmse = np.sqrt(np.sum((self.im2 - self.im1) ** 2) / float(self.im1.size))
        max_val = max(np.max(self.im1), np.max(self.im2))
        min_val = min(np.min(self.im1), np.min(self.im2))
        return 1 - (rmse / (max_val - min_val))
    def rmse(self):
        """ Root Mean Square Error """
        rmse = np.sqrt(np.sum((self.im1 - self.im2) ** 2) / float(self.im1.size))
        return rmse
    def ssim(self, window, k=(0.01, 0.03), l=255):
        from scipy.signal import fftconvolve
        """See https://ece.uwaterloo.ca/~z70wang/research/ssim/"""
        # Check if the window is smaller than the images.
        for a, b in zip(window.shape, self.im1.shape):
            if a > b:
                return None, None
        # Values in k must be positive according to the base implementation.
        for ki in k:
            if ki < 0:
                return None, None
    
        c1 = (k[0] * l) ** 2
        c2 = (k[1] * l) ** 2
        window = window/np.sum(window)
    
        mu1 = fftconvolve(self.im1, window, mode='valid')
        mu2 = fftconvolve(self.im2, window, mode='valid')
        mu1_sq = mu1 * mu1
        mu2_sq = mu2 * mu2
        mu1_mu2 = mu1 * mu2
        sigma1_sq = fftconvolve(self.im1 * self.im1, window, mode='valid') - mu1_sq
        sigma2_sq = fftconvolve(self.im2 * self.im2, window, mode='valid') - mu2_sq
        sigma12 = fftconvolve(self.im1 * self.im2, window, mode='valid') - mu1_mu2
    
        if c1 > 0 and c2 > 0:
            num = (2 * mu1_mu2 + c1) * (2 * sigma12 + c2)
            den = (mu1_sq + mu2_sq + c1) * (sigma1_sq + sigma2_sq + c2)
            ssim_map = num / den
        else:
            num1 = 2 * mu1_mu2 + c1
            num2 = 2 * sigma12 + c2
            den1 = mu1_sq + mu2_sq + c1
            den2 = sigma1_sq + sigma2_sq + c2
            ssim_map = np.ones(np.shape(mu1))
            index = (den1 * den2) > 0
            ssim_map[index] = (num1[index] * num2[index]) / (den1[index] * den2[index])
            index = (den1 != 0) & (den2 == 0)
            ssim_map[index] = num1[index] / den1[index]    
        mssim = ssim_map.mean()
        return mssim, ssim_map

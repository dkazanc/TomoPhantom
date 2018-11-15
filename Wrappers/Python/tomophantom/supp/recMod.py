#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
A reconstruction class for TomoPhantom, currently includes direct methods:
-- Fourier Slice Theorem reconstruction (adopted from Tim Day's code)
-- Filtered Back Projection (to be added, wip)

@author: Daniil Kazantsev
"""

import numpy as np

class RecTools:
    """ Class for reconstruction using TomoPhantom """
    def __init__(self, DetectorsDim, AnglesVec, ObjSize):
        self.DetectorsDim = DetectorsDim # detector dimension
        self.AnglesVec = AnglesVec # angles arrray in radians
        if ObjSize is tuple: 
            raise (" Reconstruction is currently for square or cubic objects only ")
        else:
            self.ObjSize = ObjSize # size of the object
        
    def fourier(self, sinogram, method='linear'):
        """ 
        2D Reconstruction using Fourier slice theorem (scipy required) 
        for griddata interpolation module choose nearest, linear or cubic
        """
        if sinogram.ndim == 3:
            raise ("Fourier method is currently for 2D data only, use FBP if 3D needed ")
        else:
            pass
        if ((method == 'linear') or (method == 'nearest') or (method == 'cubic')):
            pass
        else:
            raise ("For griddata interpolation module choose nearest, linear or cubic ")
        import scipy.interpolate
        import scipy.fftpack
        import scipy.misc
        import scipy.ndimage.interpolation
        
        # Fourier transform the rows of the sinogram, move the DC component to the row's centre
        sinogram_fft_rows=scipy.fftpack.fftshift(scipy.fftpack.fft(scipy.fftpack.ifftshift(sinogram,axes=1)),axes=1)
        
        """
        V  = 100
        plt.figure()
        plt.subplot(121)
        plt.title("Sinogram rows FFT (real)")
        plt.imshow(np.real(sinogram_fft_rows),vmin=-V,vmax=V)
        plt.subplot(122)
        plt.title("Sinogram rows FFT (imag)")
        plt.imshow(np.imag(sinogram_fft_rows),vmin=-V,vmax=V)
        """
        # Coordinates of sinogram FFT-ed rows' samples in 2D FFT space
        a = -self.AnglesVec
        r=np.arange(self.DetectorsDim) - self.DetectorsDim/2
        r,a=np.meshgrid(r,a)
        r=r.flatten()
        a=a.flatten()
        srcx=(self.DetectorsDim /2)+r*np.cos(a)
        srcy=(self.DetectorsDim /2)+r*np.sin(a)
        
        # Coordinates of regular grid in 2D FFT space
        dstx,dsty=np.meshgrid(np.arange(self.DetectorsDim),np.arange(self.DetectorsDim))
        dstx=dstx.flatten()
        dsty=dsty.flatten()
        
        """
        V = 100
        plt.figure()
        plt.title("Sinogram samples in 2D FFT (abs)")
        plt.scatter(srcx, srcy,c=np.absolute(sinogram_fft_rows.flatten()), marker='.', edgecolor='none', vmin=-V, vmax=V)
        """
        # Interpolate the 2D Fourier space grid from the transformed sinogram rows
        fft2=scipy.interpolate.griddata((srcy,srcx), sinogram_fft_rows.flatten(), (dsty,dstx), method, fill_value=0.0).reshape((self.DetectorsDim,self.DetectorsDim))
        """
        plt.figure()
        plt.suptitle("FFT2 space")
        plt.subplot(221)
        plt.title("Recon (real)")
        plt.imshow(np.real(fft2),vmin=-V,vmax=V)
        plt.subplot(222)
        plt.title("Recon (imag)")
        plt.imshow(np.imag(fft2),vmin=-V,vmax=V)
        """
        
        """
        # Show 2D FFT of target, just for comparison
        expected_fft2=scipy.fftpack.fftshift(scipy.fftpack.fft2(scipy.fftpack.ifftshift(phantom_2D)))
        
        plt.subplot(223)
        plt.title("Expected (real)")
        plt.imshow(np.real(expected_fft2),vmin=-V,vmax=V)
        plt.subplot(224)
        plt.title("Expected (imag)")
        plt.imshow(np.imag(expected_fft2),vmin=-V,vmax=V)
        """
        # Transform from 2D Fourier space back to a reconstruction of the target
        recon=np.real(scipy.fftpack.fftshift(scipy.fftpack.ifft2(scipy.fftpack.ifftshift(fft2))))
        
        # Cropping reconstruction to size of the original image
        image = recon[int(((self.DetectorsDim-self.ObjSize)/2)+1):self.DetectorsDim-int(((self.DetectorsDim-self.ObjSize)/2)-1),int(((self.DetectorsDim-self.ObjSize)/2)):self.DetectorsDim-int(((self.DetectorsDim-self.ObjSize)/2))]
        return image
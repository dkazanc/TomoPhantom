#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 12:28:29 2018
A reconstruction class for TomoPhantom, currently includes methods:
-- Fourier Slice Theorem reconstruction (adopted from Tim Day's code)
-- FISTA: FISTA algorithm with regularisation or not
-- Filtered Back Projection (wip)

@author: Daniil Kazantsev
"""
def powermethod(operator):
    import numpy as np
    from numpy import linalg as LA
    niter = 15
    x1 = np.float32(np.random.randn(operator.ObjSize,operator.ObjSize))
    y = operator.forwproj(x1)
    for iter in range(0,niter):
        x1 = operator.backproj(y)
        s = LA.norm(x1)
        x1 = x1/s
        y = operator.forwproj(x1)
    return s

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
    def FISTA(self, sinogram, datafidelity='LS', regularisation='none', iterationsFISTA=100, tolerance = 1e-06, device='gpu'):
        """ 
        Reconstruction method based on A. Beck and M. Teboulle, 
        A fast iterative shrinkage-thresholding algorithm for linear inverse problems,
        SIAM Journal on Imaging Sciences, vol. 2, no. 1, pp. 183–202, 2009.
        
        If regularisation is used install: CCPi-RGL toolkit 
        conda install ccpi-regulariser -c ccpi -c conda-forge
        https://github.com/vais-ral/CCPi-Regularisation-Toolkit
        """
        if sinogram.ndim == 2:
            # 2D reconstruction
            from tomophantom.supp.astraOP import AstraTools
            X = np.zeros((self.ObjSize,self.ObjSize), 'float32')
            
            Atools = AstraTools(self.DetectorsDim, self.AnglesVec, self.ObjSize, device) # initiate ASTRA class object
            X = np.zeros((self.ObjSize,self.ObjSize), 'float32')
            t = 1.0
            X_t = np.copy(X)
            L_const_inv = 1.0/powermethod(Atools) # get Lipschitz constant
            
            # Outer FISTA iterations
            for iter in range(0,iterationsFISTA):
                X_old = X
                t_old = t
                if (datafidelity == 'LS'):
                    grad_fidelity = Atools.backproj(Atools.forwproj(X_t) - sinogram) # gradient step for the LS fidelity
                else:
                    raise ("Choose data fidelity term as 'LS', 'PWLS'")
                X = X_t - L_const_inv*grad_fidelity
                # <<<<<< ADD HERE >>>>>>> the proximal operator of the regulariser
                t = (1.0 + np.sqrt(1.0 + 4.0*t**2))*0.5; # updating t
                X_t = X + ((t_old - 1.0)/t)*(X - X_old) # updating X
        if sinogram.ndim == 3:
            raise ("TODO")
        return X
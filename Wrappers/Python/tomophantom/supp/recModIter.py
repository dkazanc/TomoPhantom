#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
A reconstruction class for TomoPhantom, currently includes following iterative methods:
-- Regularised FISTA algorithm

@author: Daniil Kazantsev
"""
def powermethod(operator):
    # power iteration algorithm to  calculate the eigenvalue of the operator (projection matrix)
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
from numpy import linalg as LA

class RecTools:
    """ 
    A class for iterative reconstruction algorithms using ASTRA (only for forward/back operations)
    If regularisation is used install: CCPi-RGL toolkit 
    conda install ccpi-regulariser -c ccpi -c conda-forge
    https://github.com/vais-ral/CCPi-Regularisation-Toolkit
    """
    def __init__(self, 
              DetectorsDimH,  # DetectorsDimH # detector dimension (horizontal)
              DetectorsDimV,  # DetectorsDimV # detector dimension (vertical) for 3D case only
              AnglesVec, # array of angles in radians
              ObjSize, # a scalar to define reconstructed object dimensions
              IterativeMethod, # iterative method
              datafidelity,# data fidelity, choose LS, PWLS (wip), GH (wip), Student (wip)
              iterationsOuter, # the number of OUTER iterations
              tolerance, # tolerance to stop outer iterations earlier
              device):
        if ObjSize is tuple: 
            raise (" Reconstruction is currently available for square or cubic objects only ")
        else:
            self.ObjSize = ObjSize # size of the object
        
        self.tolerance = tolerance
        self.iterationsOuter = iterationsOuter
        self.datafidelity = datafidelity
        if device is None:
            self.device = 'gpu'
        else:
            self.device = device
        
        if DetectorsDimV is None:
            # Creating Astra class specific to 2D parallel geometry
            from tomophantom.supp.astraOP import AstraTools
            self.Atools = AstraTools(DetectorsDimH, AnglesVec, ObjSize, device) # initiate 2D ASTRA class object
            self.geom = '2D'
        else:
            # Creating Astra class specific to 3D parallel geometry
            from tomophantom.supp.astraOP import AstraTools3D
            self.Atools = AstraTools3D(DetectorsDimH, DetectorsDimV, AnglesVec, ObjSize) # initiate 3D ASTRA class object
            self.geom = '3D'

        if (IterativeMethod == 'FISTA'):
            # set all parameters here for FISTA
            self.L_const_inv = 1.0/powermethod(self.Atools) # get Lipschitz constant

    def FISTA(self, projdata, InitialObject = 0, regularisation = None, regularisation_parameter = 0.01, regularisation_iterations = 50):
        """ 
        Reconstruction method based on A. Beck and M. Teboulle, 
        A fast iterative shrinkage-thresholding algorithm for linear inverse problems,
        SIAM Journal on Imaging Sciences, vol. 2, no. 1, pp. 183–202, 2009.
        
        Input parameters:
        projdata - tomographic projection data in 2D (sinogram) or 3D array
        InitialObject - initialise reconstruction with an array
        regularisation, # enable regularisation  with CCPi - RGL toolkit
        regularisation_parameter, # regularisation parameter if regularisation is not None
        regularisation_iterations, # the number of INNER iterations for regularisation
        
        """
        if (self.geom == '2D'):
            # 2D reconstruction
            # initialise the solution
            if (np.size(InitialObject) == self.ObjSize**2):
                # the object has been initialised with an array
                X = InitialObject
                del InitialObject
            else:
                X = np.zeros((self.ObjSize,self.ObjSize), 'float32')
            # If the regularisation is used, then the dependency is on the CCPi-RGL toolkit
            if (regularisation == 'ROF_TV'):
                # Rudin - Osher - Fatemi Total variation method
                from ccpi.filters.regularisers import ROF_TV
                # setting some default values
                regularisation_iterations = 250
                time_marching_parameter = 0.0025 # gradient step parameter
            if (regularisation == 'FGP_TV'):
                # Fast-Gradient-Projection Total variation method
                from ccpi.filters.regularisers import FGP_TV
                # setting some default values
                regularisation_iterations = 100
                tolerance_regul = 1e-06 # tolerance to stop regularisation
                methodTV = 0 # 0/1 - isotropic/anisotropic TV
                nonneg = 0 # 0/1 disabled/enabled nonnegativity
            if (regularisation == 'SB_TV'):
                # Split Bregman Total variation method
                from ccpi.filters.regularisers import SB_TV
                # setting some default values
                regularisation_iterations = 50
                tolerance_regul = 1e-06 # tolerance to stop regularisation
                methodTV = 0 # 0/1 - isotropic/anisotropic TV
                
        if (self.geom == '3D'):
            # initialise the solution
            if (np.size(InitialObject) == self.ObjSize**3):
                # the object has been initialised with an array
                X = InitialObject
                del InitialObject
            else:
                X = np.zeros((self.ObjSize,self.ObjSize,self.ObjSize), 'float32')
#****************************************************************************#
            # FISTA algorithm begins here:
            t = 1.0
            denomN = 1.0/np.size(X)
            X_t = np.copy(X)
            
            # Outer FISTA iterations
            for iter in range(0,self.iterationsOuter):
                X_old = X
                t_old = t
                if (self.datafidelity == 'LS'):
                    grad_fidelity = self.Atools.backproj(self.Atools.forwproj(X_t) - projdata) # gradient step for the LS fidelity
                else:
                    raise ("Choose the data fidelity term as 'LS', 'PWLS'")
                X = X_t - self.L_const_inv*grad_fidelity
                # stopping criteria
                nrm = LA.norm(X - X_old)*denomN
                if nrm > self.tolerance:
                    # The proximal operator of the chosen regulariser
                    if (regularisation == 'ROF_TV'):
                        X = ROF_TV(X, regularisation_parameter, regularisation_iterations, time_marching_parameter, self.device)
                    if (regularisation == 'FGP_TV'):
                        X = FGP_TV(X, regularisation_parameter, regularisation_iterations, tolerance_regul, methodTV, nonneg, 0, self.device)
                    if (regularisation == 'SB_TV'):
                        X = SB_TV(X, regularisation_parameter, regularisation_iterations, tolerance_regul, methodTV, 0, self.device)
                    t = (1.0 + np.sqrt(1.0 + 4.0*t**2))*0.5; # updating t variable
                    X_t = X + ((t_old - 1.0)/t)*(X - X_old) # updating X
#****************************************************************************#
        return X
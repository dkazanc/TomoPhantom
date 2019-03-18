#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D analytical phantoms and their sinograms with added noise and artifacts
Sinograms then reconstructed using ASTRA TOOLBOX 

>>>>> Dependencies (reconstruction): <<<<<
1. ASTRA toolbox: conda install -c astra-toolbox astra-toolbox
2. TomoRec: conda install -c dkazanc tomorec
or install from https://github.com/dkazanc/TomoRec

This demo demonstrates frequent inaccuracies which are accosiated with X-ray imaging:
zingers, rings and noise

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP2D
import os
import tomophantom
from tomophantom.supp.qualitymetrics import QualityTools

model = 4 # select a model
N_size = 512 # set dimension of the phantom
# one can specify an exact path to the parameters file
# path_library2D = '../../../PhantomLibrary/models/Phantom2DLibrary.dat'
path = os.path.dirname(tomophantom.__file__)
path_library2D = os.path.join(path, "Phantom2DLibrary.dat")
phantom_2D = TomoP2D.Model(model, N_size, path_library2D)

plt.close('all')
plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.imshow(phantom_2D, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}''{}'.format('2D Phantom using model no.',model))

# create sinogram analytically
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0.0,179.9,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180.0)
P = int(np.sqrt(2)*N_size) #detectors

sino_an = TomoP2D.ModelSino(model, N_size, P, angles, path_library2D)

plt.figure(2)
plt.rcParams.update({'font.size': 21})
plt.imshow(sino_an,  cmap="gray")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Analytical sinogram of model no.',model))

#%%
# Adding artifacts and noise
from tomophantom.supp.artifacts import ArtifactsClass

# adding noise
artifacts_add = ArtifactsClass(sino_an)
#noisy_sino = artifacts_add.noise(sigma=0.1,noisetype='Gaussian')
noisy_sino = artifacts_add.noise(sigma=10000,noisetype='Poisson')

# adding object shifts (misalignment)
artifacts_add = ArtifactsClass(noisy_sino)
noisy_sino_misalign = artifacts_add.shifts(maxamplitude = 10)

# adding zingers
artifacts_add =ArtifactsClass(noisy_sino)
noisy_zing = artifacts_add.zingers(percentage=0.25, modulus = 10)

#adding stripes
artifacts_add =ArtifactsClass(noisy_zing)
noisy_zing_stripe = artifacts_add.stripes(percentage=1, maxthickness = 1)
noisy_zing_stripe[noisy_zing_stripe < 0] = 0

plt.figure()
plt.rcParams.update({'font.size': 21})
plt.imshow(noisy_zing_stripe,cmap="gray")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Analytical noisy sinogram with artifacts.',model))
#%%
# initialise TomoRec DIRECT reconstruction class ONCE
from tomorec.methodsDIR import RecToolsDIR
RectoolsDIR = RecToolsDIR(DetectorsDimH = P,  # DetectorsDimH # detector dimension (horizontal)
                    DetectorsDimV = None,  # DetectorsDimV # detector dimension (vertical) for 3D case only
                    AnglesVec = angles_rad, # array of angles in radians
                    ObjSize = N_size, # a scalar to define reconstructed object dimensions
                    device='cpu')
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using Fourier Slice method")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
RecFourier = RectoolsDIR.fourier(sino_an,'linear') 
plt.figure() 
plt.imshow(RecFourier, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Fourier slice reconstruction')
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using FBP (TomoRec)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
FBPrec_ideal = RectoolsDIR.FBP(sino_an)  # ideal reconstruction
FBPrec_error = RectoolsDIR.FBP(noisy_zing_stripe) # reconstruction with artifacts
FBPrec_misalign = RectoolsDIR.FBP(noisy_sino_misalign) # reconstruction with misalignment

plt.figure()
plt.subplot(131)
plt.imshow(FBPrec_ideal, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Ideal FBP reconstruction')
plt.subplot(132)
plt.imshow(FBPrec_error, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Erroneous data FBP Reconstruction')
plt.subplot(133)
plt.imshow(FBPrec_misalign, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Misaligned noisy FBP Reconstruction')
plt.show()

plt.figure() 
plt.imshow(abs(FBPrec_ideal-FBPrec_error), vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP reconsrtuction differences')
#%%
# initialise TomoRec ITERATIVE reconstruction class ONCE
from tomorec.methodsIR import RecToolsIR
RectoolsIR = RecToolsIR(DetectorsDimH = P,  # DetectorsDimH # detector dimension (horizontal)
                    DetectorsDimV = None,  # DetectorsDimV # detector dimension (vertical) for 3D case only
                    AnglesVec = angles_rad, # array of angles in radians
                    ObjSize = N_size, # a scalar to define reconstructed object dimensions
                    datafidelity='LS',# data fidelity, choose LS, PWLS (wip), GH (wip), Student (wip)
                    nonnegativity='ENABLE', # enable nonnegativity constraint (set to 'ENABLE')
                    OS_number = None, # the number of subsets, NONE/(or > 1) ~ classical / ordered subsets
                    tolerance = 1e-06, # tolerance to stop outer iterations earlier
                    device='gpu')
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using SIRT (TomoRec)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
iterationsSIRT = 250
SIRTrec_ideal = RectoolsIR.SIRT(sino_an,iterationsSIRT) # ideal reconstruction
SIRTrec_error = RectoolsIR.SIRT(noisy_zing_stripe,iterationsSIRT) # error reconstruction

plt.figure()
plt.subplot(121)
plt.imshow(SIRTrec_ideal, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Ideal SIRT reconstruction (ASTRA)')
plt.subplot(122)
plt.imshow(SIRTrec_error, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Erroneous data SIRT Reconstruction (ASTRA)')
plt.show()

plt.figure() 
plt.imshow(abs(SIRTrec_ideal-SIRTrec_error), vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('SIRT reconsrtuction differences')
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing using FISTA method (TomoRec)")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
lc = RectoolsIR.powermethod() # calculate Lipschitz constant

# Run FISTA reconstrucion algorithm without regularisation
RecFISTA = RectoolsIR.FISTA(noisy_zing_stripe, iterationsFISTA = 150, lipschitz_const = lc)

# Run FISTA reconstrucion algorithm with regularisation 
RecFISTA_reg = RectoolsIR.FISTA(noisy_zing_stripe, iterationsFISTA = 150, regularisation = 'ROF_TV', regularisation_parameter = 0.003, lipschitz_const = lc)

plt.figure()
plt.subplot(121)
plt.imshow(RecFISTA, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FISTA reconstruction')
plt.subplot(122)
plt.imshow(RecFISTA_reg, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('TV-Regularised FISTA reconstruction')
plt.show()

# calculate errors 
Qtools = QualityTools(phantom_2D, RecFISTA)
RMSE_FISTA = Qtools.rmse()
Qtools = QualityTools(phantom_2D, RecFISTA_reg)
RMSE_FISTA_reg = Qtools.rmse()
print("RMSE for FISTA is {}".format(RMSE_FISTA))
print("RMSE for regularised FISTA is {}".format(RMSE_FISTA_reg))
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing using ADMM method (TomoRec)")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
# Run ADMM reconstrucion algorithm with regularisation 
RecADMM_reg = RectoolsIR.ADMM(noisy_zing_stripe, rho_const = 5000.0, iterationsADMM = 10, regularisation = 'ROF_TV', regularisation_parameter = 0.3)

plt.figure()
plt.imshow(RecADMM_reg, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('ADMM reconstruction')
plt.show()

# calculate errors 
Qtools = QualityTools(phantom_2D, RecADMM_reg)
RMSE_ADMM_reg = Qtools.rmse()
print("RMSE for regularised ADMM is {}".format(RMSE_ADMM_reg))
#%%

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D analytical phantoms and their sinograms with added noise and artifacts
Sinograms then reconstructed using ASTRA TOOLBOX and TomoPy packages

>>>> Prerequisites: ASTRA toolbox and TomoPy packages <<<<<

This demo demonstrates frequent inaccuracies which are accosiated with X-ray imaging:
zingers, rings and noise

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP2D

model = 4 # select a model
N_size = 512
#specify a full path to the parameters file
pathTP = '../../functions/models/Phantom2DLibrary.dat'
#objlist = modelfile2Dtolist(pathTP, model) # one can extract parameters
#This will generate a N_size x N_size phantom (2D)
phantom_2D = TomoP2D.Model(model, N_size, pathTP)

plt.close('all')
plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.imshow(phantom_2D, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}''{}'.format('2D Phantom using model no.',model))

# create sinogram analytically
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors

sino_an = TomoP2D.ModelSino(model, N_size, P, angles, pathTP)

plt.figure(2)
plt.rcParams.update({'font.size': 21})
plt.imshow(sino_an,  cmap="gray")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Analytical sinogram of model no.',model))

#%%
# Adding artifacts and noise
from tomophantom.supp.artifacts import ArtifactsClass

# adding noise
artifacts_add =ArtifactsClass(sino_an)
#noisy_sino = artifacts_add.noise(sigma=0.1,noisetype='Gaussian')
noisy_sino = artifacts_add.noise(sigma=10000,noisetype='Poisson')

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
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using FBP (ASTRA-TOOLBOX)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
from tomophantom.supp.astraOP import AstraTools

Atools = AstraTools(P, angles_rad, N_size, 'cpu') # initiate a class object

FBPrec_ideal = Atools.fbp2D(sino_an) # ideal reconstruction
FBPrec_error = Atools.fbp2D(noisy_zing_stripe) # error reconstruction

plt.figure()
plt.subplot(121)
plt.imshow(FBPrec_ideal, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Ideal reconstruction (ASTRA)')
plt.subplot(122)
plt.imshow(FBPrec_error, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Erroneous Reconstruction (ASTRA)')
plt.show()

plt.figure() 
plt.imshow(abs(FBPrec_ideal-FBPrec_error), vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP reconsrtuction differences')
rmse2 = np.linalg.norm(FBPrec_ideal-FBPrec_error)/np.linalg.norm(FBPrec_error)

#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using gridrec (TomoPy)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
import tomopy

sinoTP = np.zeros((angles_num,1,P),dtype='float32')
sinoTP[:,0,:] = sino_an
rot_center = tomopy.find_center(sinoTP, angles_rad, init=290, ind=0, tol=0.5)
reconTomoPy_ideal = tomopy.recon(sinoTP, angles_rad, center=rot_center, algorithm='gridrec')
sinoTP[:,0,:] = noisy_zing_stripe
reconTomoPy_noisy = tomopy.recon(sinoTP, angles_rad, center=rot_center, algorithm='gridrec')


plt.figure()
plt.subplot(121)
plt.imshow(reconTomoPy_ideal[0,:,:], vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Ideal reconstruction (TomoPy)')
plt.subplot(122)
plt.imshow(reconTomoPy_noisy[0,:,:], vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Erroneous Reconstruction (TomoPy)')
plt.show()
#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D analytical phantoms and their sinograms
If one needs to modify/add phantoms, please edit Phantom2DLibrary.dat
>>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

!Run script from "Demos" folder in order to ensure a correct path to *dat file!

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP2D

model = 1 # selecting a model
N_size = 512
#specify a full path to the parameters file
pathTP = '../../functions/models/Phantom2DLibrary.dat'
#objlist = modelfile2Dtolist(pathTP, model) # one can extract parameters
#This will generate a N_size x N_size phantom (2D)
phantom_2D = TomoP2D.Model(model, N_size, pathTP)

plt.close('all')
plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.imshow(phantom_2D, vmin=0, vmax=1, cmap="BuPu")
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
plt.imshow(sino_an,  cmap="BuPu")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Analytical sinogram of model no.',model))
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using FBP (ASTRA-TOOLBOX)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
from tomophantom.supp.astraOP import AstraTools

Atools = AstraTools(P, angles_rad, N_size, 'cpu') # initiate a class object
sino_num_ASTRA = Atools.forwproj(phantom_2D) # generate numerical sino (Ax)
#x = Atools.backproj(sino_an) # generate backprojection (A'b)

plt.figure(2) 
plt.subplot(121)
plt.imshow(sino_an,cmap="BuPu")
plt.title('Analytical sinogram')
plt.subplot(122)
plt.imshow(sino_num_ASTRA,cmap="BuPu")
plt.title('Numerical sinogram')
plt.show()
#calculate norm
rmse1 = np.linalg.norm(sino_an - sino_num_ASTRA)/np.linalg.norm(sino_num_ASTRA)

print ("Reconstructing analytical sinogram using FBP (astra TB)...")
FBPrec1 = Atools.fbp2D(sino_an)
plt.figure(3) 
plt.imshow(FBPrec1, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP Reconstructed Phantom (analyt)')

print ("Reconstructing numerical sinogram using FBP (astra TB)...")
FBPrec2 = Atools.fbp2D(sino_num_ASTRA)

plt.figure(4) 
plt.imshow(FBPrec2, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP Reconstructed Phantom (numeric)')

plt.figure(5) 
plt.imshow(abs(FBPrec1-FBPrec2), vmin=0, vmax=0.05, cmap="BuPu")
plt.colorbar(ticks=[0, 0.02, 0.05], orientation='vertical')
plt.title('FBP rec differences')
rmse2 = np.linalg.norm(FBPrec1 - FBPrec2)/np.linalg.norm(FBPrec2)
#%%

#%%

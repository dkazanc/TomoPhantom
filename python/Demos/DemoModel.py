#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D/3D analytical phantoms and their sinograms
If one needs to modify/add phantoms, please edit Phantom2DLibrary.dat or
Phantom3DLibrary.dat
>>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

Run demo from the folder "Demos"

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP2D
from astraOP import AstraTools
#%%
model = 4
N_size = 512
#specify a full path to the parameters file
pathTP = '../../functions/models/Phantom2DLibrary.dat'
#This will generate a N_size x N_size phantom (2D)
phantom_2D = TomoP2D.Model(model, N_size, pathTP)

plt.close('all')
plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.imshow(phantom_2D, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}''{}'.format('2D Phantom using model no.',model))
#%%
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
Atools = AstraTools(P, angles_rad - 0.5*np.pi, N_size, 'cpu') # initiate a class object
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
#%%
print ("Reconstructing analytical sinogram using FBP (astra TB)...")
FBPrec1 = Atools.fbp2D(sino_an)

plt.figure(3) 
plt.imshow(FBPrec1, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP Reconstructed Phantom (analyt)')
#%%
print ("Reconstructing numerical sinogram using FBP (astra TB)...")
FBPrec2 = Atools.fbp2D(sino_num_ASTRA)

plt.figure(4) 
plt.imshow(FBPrec2, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP Reconstructed Phantom (numeric)')
#%%
plt.figure(5) 
plt.imshow(abs(FBPrec1-FBPrec2), vmin=0, vmax=0.05, cmap="BuPu")
plt.colorbar(ticks=[0, 0.02, 0.05], orientation='vertical')
plt.title('FBP rec differences')
#%%
"""
print ("Reconstructing using SIRT...")
SIRTrec = Atools.sirt2D(sino_an, 100)

plt.figure(4) 
plt.imshow(SIRTrec, vmin=0, vmax=1,cmap="BuPu")
plt.title('SIRT Reconstructed Phantom')
"""
#%%
"""
import timeit
from tomophantom import TomoP3D
import matplotlib.pyplot as plt

print ("Building 3D phantom using TomoPhantom software")
tic=timeit.default_timer()
model = 1
N_size = 256
#specify a full path to the parameters file
pathTP3 = '../../functions/models/Phantom3DLibrary.dat'
#This will generate a N_size x N_size x N_size phantom (3D)
phantom_tm = TomoP3D.Model(model, N_size, pathTP3)
toc=timeit.default_timer()
Run_time = toc - tic
print("Phantom has been built in {} seconds".format(Run_time))

sliceSel = int(0.5*N_size)
#plt.gray()
plt.figure(5) 
plt.subplot(131)
plt.imshow(phantom_tm[sliceSel,:,:],vmin=0, vmax=1)
plt.title('3D Phantom, axial view')

plt.subplot(132)
plt.imshow(phantom_tm[:,sliceSel,:],vmin=0, vmax=1)
plt.title('3D Phantom, coronal view')

plt.subplot(133)
plt.imshow(phantom_tm[:,:,sliceSel],vmin=0, vmax=1)
plt.title('3D Phantom, sagittal view')
plt.show()
"""

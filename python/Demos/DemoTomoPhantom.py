#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D/3D analytical phantoms and their sinograms
If one needs to modify/add phantoms, please edit Phantom2DLibrary.dat or
Phantom3DLibrary.dat
>>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import phantom2d
from tomophantom import phantom3d
from astraOP import AstraTools
#%%
model = 11
N_size = 512
#specify a full path to the parameters file
pathTP = '/home/algol/Documents/DEV/TomoPhantom/functions/models/Phantom2DLibrary.dat'
#This will generate a N_size x N_size phantom (2D)
phantom_2D = phantom2d.buildPhantom2D(model, N_size, pathTP)

plt.figure(1)
plt.rcParams.update({'font.size': 22})
plt.imshow(phantom_2D, vmin=0, vmax=1)
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}''{}'.format('2D Phantom using model no.',model))
#%%
# create sinogram analytically
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors

sino_an = phantom2d.buildSino2D(model, N_size, P, angles, pathTP, 1)

plt.figure(2)
plt.rcParams.update({'font.size': 22})
plt.imshow(sino_an)
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Analytical sinogram of model no.',model))
#%%
Atools = AstraTools(P, angles_rad - 0.5*np.pi, N_size, 'cpu') # initiate a class object
sino_num_ASTRA = Atools.forwproj(phantom_2D) # generate numerical sino (Ax)
#x = Atools.backproj(sino_an) # generate backprojection (A'b)

plt.figure(2) 
plt.subplot(121)
plt.imshow(sino_an)
plt.title('Analytical sinogram')
plt.subplot(122)
plt.imshow(sino_num_ASTRA)
plt.title('Numerical sinogram')
plt.show()
#%%
print ("Reconstructing using FBP...")
FBPrec = Atools.fbp2D(sino_an)

plt.figure(3) 
plt.imshow(FBPrec, vmin=0, vmax=1)
plt.title('FBP Reconstructed Phantom')
#%%
print ("Reconstructing using SIRT...")
SIRTrec = Atools.sirt2D(sino_an, 100)

plt.figure(4) 
plt.imshow(SIRTrec, vmin=0, vmax=1)
plt.title('SIRT Reconstructed Phantom')

#%%
import timeit
print ("Building 3D phantom using TomoPhantom software")
tic=timeit.default_timer()
model = 9
N_size = 512
#specify a full path to the parameters file
pathTP3 = '/home/algol/Documents/DEV/TomoPhantom/functions/models/Phantom3DLibrary.dat'
#This will generate a N_size x N_size x N_size phantom (3D)
phantom_tm = phantom3d.buildPhantom3D(model, N_size, pathTP3)
toc=timeit.default_timer()
Run_time = toc - tic
print("Phantom has been built in {} seconds".format(Run_time))

sliceSel = int(0.5*N_size)
#plt.gray()
plt.figure(5) 
plt.subplot(121)
plt.imshow(phantom_tm[sliceSel,:,:],vmin=0, vmax=1)
plt.title('3D Phantom, axial view')

plt.subplot(122)
plt.imshow(phantom_tm[:,sliceSel,:],vmin=0, vmax=1)
plt.title('3D Phantom, coronal view')
plt.show()


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D/3D analytical phantoms and their sinograms
If one needs to modify/add phantoms, please edit Phantom2DLibrary.dat or
Phantom3DLibrary.dat
>>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

Main difference from DemoModel.py is that we extract all parameters from the 
library file using Python and then pass it to the Object function instead of model.
This can be helpful if one would like to avoid using library files and can
pass parameters directly into object function

!Run script from "Demos" folder in order to ensure a correct path to *dat file!

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP2D
from libraryToDict import modelfile2Dtolist

model = 11
N_size = 512
#specify a full path to the parameters file
pathTP = '../../functions/models/Phantom2DLibrary.dat'
objlist = modelfile2Dtolist(pathTP, model) # extract parameters using Python
#This will generate a N_size x N_size phantom (2D)
phantom_2D = TomoP2D.Object(N_size, objlist)

plt.close('all')
plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.imshow(phantom_2D, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}''{}'.format('2D Phantom using model no.',model))

# create sinogram analytically
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0.0,179.9,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180.0)
P = int(np.sqrt(2)*N_size) #detectors

sino_an = TomoP2D.ObjectSino(N_size, P, angles, objlist)

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

plt.figure(3) 
plt.subplot(121)
plt.imshow(sino_an,cmap="BuPu")
plt.title('Analytical sinogram')
plt.subplot(122)
plt.imshow(sino_num_ASTRA,cmap="BuPu")
plt.title('Numerical sinogram')
plt.show()

print ("Reconstructing analytical sinogram using FBP (astra TB)...")
FBPrec1 = Atools.fbp2D(sino_an)
plt.figure(4) 
plt.imshow(FBPrec1, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP Reconstructed Phantom (analyt)')

print ("Reconstructing numerical sinogram using FBP (astra TB)...")
FBPrec2 = Atools.fbp2D(sino_num_ASTRA)
plt.figure(5) 
plt.imshow(FBPrec2, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP Reconstructed Phantom (numeric)')

plt.figure(6) 
plt.imshow(abs(FBPrec1-FBPrec2), vmin=0, vmax=0.05, cmap="BuPu")
plt.colorbar(ticks=[0, 0.02, 0.05], orientation='vertical')
plt.title('FBP rec differences')
rmse2 = np.linalg.norm(FBPrec1 - FBPrec2)/np.linalg.norm(FBPrec2)
#%%
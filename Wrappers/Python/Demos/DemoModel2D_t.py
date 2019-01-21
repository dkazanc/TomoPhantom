#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate temporal (2D + time)  analytical phantoms and their sinograms
If one needs to modify/add phantoms, please edit Phantom2DLibrary.dat or
Phantom3DLibrary.dat
! Note that all temporal phantoms start from no. 100 

>>>>> Optional dependencies (reconstruction mainly): <<<<<
1. ASTRA toolbox: conda install -c astra-toolbox astra-toolbox
2. TomoRec: conda install -c dkazanc tomorec
or install from https://github.com/dkazanc/TomoRec


@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import tomophantom
from tomophantom import TomoP2D

model = 102  # note that the selected model is temporal (2D + time)
N_size = 512 # set dimension of the phantom
# one can specify an exact path to the parameters file
# path_library2D = '../../../PhantomLibrary/models/Phantom2DLibrary.dat'
path = os.path.dirname(tomophantom.__file__)
path_library2D = os.path.join(path, "Phantom2DLibrary.dat")
#This will generate a N_size x N_size x Time frames phantom (2D + time)
phantom_2Dt = TomoP2D.ModelTemporal(model, N_size, path_library2D)

plt.close('all')
plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.title('{}''{}'.format('2D+t phantom using model no.',model))
for sl in range(0,np.shape(phantom_2Dt)[0]):
    im = phantom_2Dt[sl,:,:]
    plt.imshow(im, vmin=0, vmax=1)
    plt.pause(.1)
    plt.draw

# create sinogram analytically
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors

sino = TomoP2D.ModelSinoTemporal(model, N_size, P, angles, path_library2D)

plt.figure(2)
plt.rcParams.update({'font.size': 21})
plt.title('{}''{}'.format('2D+t sinogram of model no.',model))
for sl in range(0,np.shape(phantom_2Dt)[0]):
    im = sino[sl,:,:].transpose()
    plt.imshow(im, vmin=0, vmax=180)
    plt.pause(.1)
    plt.draw
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using FBP (TomoRec)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
# initialise TomoRec reconstruction class ONCE
from tomorec.methodsDIR import RecToolsDIR
RectoolsDIR = RecToolsDIR(DetectorsDimH = P,  # DetectorsDimH # detector dimension (horizontal)
                    DetectorsDimV = None,  # DetectorsDimV # detector dimension (vertical) for 3D case only
                    AnglesVec = angles_rad, # array of angles in radians
                    ObjSize = N_size, # a scalar to define reconstructed object dimensions
                    device='cpu')

FBPrec = RectoolsDIR.FBP(sino[15,:,:].transpose()) # reconstruct one frame

plt.figure(3) 
plt.imshow(FBPrec, vmin=0, vmax=1)
plt.title('FBP Reconstructed Phantom')
#%%
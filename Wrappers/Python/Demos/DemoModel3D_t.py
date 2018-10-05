#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 4D (3D+time) analytical phantoms (wip: generation of 4D projection data )
If one needs to modify/add phantoms, please edit Phantom3DLibrary.dat

!Run script from "Demos" folder in order to ensure a correct path to *dat file!

@author: Daniil Kazantsev
"""
import timeit
import os
import matplotlib.pyplot as plt
import numpy as np
import tomophantom
from tomophantom import TomoP3D

print ("Building 4D phantom using TomoPhantom software")
tic=timeit.default_timer()
model = 100 # note that the selected model is temporal (3D + time)
# Define phantom dimensions using a scalar (cubic) or a tuple [N1, N2, N3]
N_size = 256 # or as a tuple of a custom size (256,256,256)
# one can specify an exact path to the parameters file
# path_library2D = '../../../PhantomLibrary/models/Phantom3DLibrary.dat'
path = os.path.dirname(tomophantom.__file__)
path_library3D = os.path.join(path, "Phantom3DLibrary.dat")
#This will generate a Time x N_size x N_size x N_size phantom (4D)
phantom_tm = TomoP3D.ModelTemporal(model, N_size, path_library3D)
toc=timeit.default_timer()
Run_time = toc - tic
print("Phantom has been built in {} seconds".format(Run_time))

for i in range(0,np.size(phantom_tm,0)): 
    sliceSel = int(0.5*N_size)
    #plt.gray()
    plt.figure(1) 
    plt.subplot(131)
    plt.imshow(phantom_tm[i,sliceSel,:,:],vmin=0, vmax=1)
    plt.title('3D Phantom, axial view')

    plt.subplot(132)
    plt.imshow(phantom_tm[i,:,sliceSel,:],vmin=0, vmax=1)
    plt.title('3D Phantom, coronal view')

    plt.subplot(133)
    plt.imshow(phantom_tm[i,:,:,sliceSel],vmin=0, vmax=1)
    plt.title('3D Phantom, sagittal view')
    plt.show()
    plt.pause(0.3)

#%%
# A capability of building a subset of vertical slices out of 4D phantom (faster)
import timeit
import os
import tomophantom
from tomophantom import TomoP3D
import matplotlib.pyplot as plt

print ("Building a subset of 3D phantom using TomoPhantom software")
tic=timeit.default_timer()
model = 101
# Define phantom dimensions using a scalar (cubic) or a tuple [Z, Y, X]
DIM = (256,256,256) # full dimension of required phantom (z, y, x)
DIM_z = (94, 158) # selected vertical subset (a slab) of the phantom
path = os.path.dirname(tomophantom.__file__)
path_library3D = os.path.join(path, "Phantom3DLibrary.dat")
#This will generate a Time x N1 x N2 x N_slab phantom (4D)
phantom_tm = TomoP3D.ModelTemporalSub(model, DIM, DIM_z, path_library3D)
#phantom_tm = TomoP3D.Model(model, DIM, pathTP3)
toc=timeit.default_timer()
Run_time = toc - tic
print("Phantom has been built in {} seconds".format(Run_time))

for i in range(0,np.size(phantom_tm,0)): 
    sliceSel = 32
    #plt.gray()
    plt.figure(1) 
    plt.subplot(131)
    plt.imshow(phantom_tm[i,sliceSel,:,:],vmin=0, vmax=1)
    plt.title('3D Phantom, axial view')

    plt.subplot(132)
    plt.imshow(phantom_tm[i,:,70,:],vmin=0, vmax=1)
    plt.title('3D Phantom, coronal view')

    plt.subplot(133)
    plt.imshow(phantom_tm[i,:,:,70],vmin=0, vmax=1)
    plt.title('3D Phantom, sagittal view')
    plt.show()
    plt.pause(0.5)
#%%
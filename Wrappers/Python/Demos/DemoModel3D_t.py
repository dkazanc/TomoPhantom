#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 4D (3D+time) analytical phantoms (wip: generation of 4D projection data )
If one needs to modify/add phantoms, please edit Phantom3DLibrary.dat
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
print ("Getting 4D projection data using TomoPhantom software")
# Projection geometry related parameters:
Horiz_det = int(np.sqrt(2)*N_size) # detector column count (horizontal)
Vert_det = N_size # detector row count (vertical) (no reason for it to be > N)
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0.0,179.9,angles_num,dtype='float32') # in degrees
angles_rad = angles*(np.pi/180.0)

projData4D_analyt= TomoP3D.ModelSinoTemporal(model, N_size, Horiz_det, Vert_det, angles, path_library3D)

time_frames = projData4D_analyt.shape[0]
intens_max = 60
sliceSel = 150
for i in range(0,time_frames):
    plt.figure(2) 
    plt.subplot(131)
    plt.imshow(projData4D_analyt[i,:,sliceSel,:],vmin=0, vmax=intens_max)
    plt.title('2D Projection (analytical)')
    plt.subplot(132)
    plt.imshow(projData4D_analyt[i,sliceSel,:,:],vmin=0, vmax=intens_max)
    plt.title('Sinogram view')
    plt.subplot(133)
    plt.imshow(projData4D_analyt[i,:,:,sliceSel],vmin=0, vmax=intens_max)
    plt.title('Tangentogram view')
    plt.show()
    plt.pause(0.3)
#%%
# A capability of building a subset of vertical slices out of 4D phantom (faster)
import timeit
import os
import tomophantom
from tomophantom import TomoP3D
import matplotlib.pyplot as plt
import numpy as np

print ("Building a subset of 4D phantom using TomoPhantom software")
tic=timeit.default_timer()
model = 101
N_size = 256 # Define phantom dimensions using a scalar value
DIM_z = (94, 158) # selected vertical subset (a slab) of the phantom
path = os.path.dirname(tomophantom.__file__)
path_library3D = os.path.join(path, "Phantom3DLibrary.dat")
phantom_tm = TomoP3D.ModelTemporalSub(model, N_size, DIM_z, path_library3D)
toc=timeit.default_timer()
Run_time = toc - tic
print("Phantom has been built in {} seconds".format(Run_time))

for i in range(0,np.size(phantom_tm,0)): 
    sliceSel = 32
    #plt.gray()
    plt.figure(1) 
    plt.subplot(131)
    plt.imshow(phantom_tm[i,sliceSel,:,:],vmin=0, vmax=1)
    plt.title('4D Phantom, axial view')

    plt.subplot(132)
    plt.imshow(phantom_tm[i,:,70,:],vmin=0, vmax=1)
    plt.title('4D Phantom, coronal view')

    plt.subplot(133)
    plt.imshow(phantom_tm[i,:,:,70],vmin=0, vmax=1)
    plt.title('4D Phantom, sagittal view')
    plt.show()
    plt.pause(0.5)
    

print ("Building a subset of 4D projection data using TomoPhantom software")
Horiz_det = int(np.sqrt(2)*N_size) # detector column count (horizontal)
Vert_det = N_size # detector row count (vertical) (no reason for it to be > N)
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0.0,179.9,angles_num,dtype='float32') # in degrees

projData4D_cut = TomoP3D.ModelSinoTemporalSub(model, N_size, Horiz_det, Vert_det, DIM_z, angles, path_library3D)

intens_max = 45
for i in range(0,np.size(projData4D_cut,0)): 
    sliceSel = 32
    #plt.gray()
    plt.figure(1) 
    plt.subplot(131)
    plt.imshow(projData4D_cut[i,sliceSel,:,:],vmin=0, vmax=intens_max)
    plt.title('Sinogram View')

    plt.subplot(132)
    plt.imshow(projData4D_cut[i,:,sliceSel,:],vmin=0, vmax=intens_max)
    plt.title('Projection view')

    plt.subplot(133)
    plt.imshow(projData4D_cut[i,:,:,sliceSel],vmin=0, vmax=intens_max)
    plt.title('Tangentogram view')
    plt.show()
    plt.pause(0.5)
#%%
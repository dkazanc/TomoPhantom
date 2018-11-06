#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

* Script to generate 3D analytical phantoms and their projection data using TomoPhantom
* Projection data can be also generated numerically and reconstructed using ASTRA TOOLBOX 

>>>>> Prerequisites: ASTRA toolbox  <<<<<
install ASTRA: conda install -c astra-toolbox astra-toolbox

@author: Daniil Kazantsev
"""
import timeit
import os
from tomophantom import TomoP3D
import matplotlib.pyplot as plt
import tomophantom

print ("Building 3D phantom using TomoPhantom software")
tic=timeit.default_timer()
model = 9 # select a model number from the library
# Define phantom dimensions using a scalar (cubic) or a tuple [N1, N2, N3]
N_size = 256 # or as a tuple of a custom size (256,256,256)
# one can specify an exact path to the parameters file
# path_library2D = '../../../PhantomLibrary/models/Phantom3DLibrary.dat'
path = os.path.dirname(tomophantom.__file__)
path_library3D = os.path.join(path, "Phantom3DLibrary.dat")
#This will generate a N_size x N_size x N_size phantom (3D) or non-cubic phantom
phantom_tm = TomoP3D.Model(model, N_size, path_library3D)
toc=timeit.default_timer()
Run_time = toc - tic
print("Phantom has been built in {} seconds".format(Run_time))

sliceSel = int(0.5*N_size)
#plt.gray()
plt.figure() 
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

import numpy as np
Horiz_det = int(np.sqrt(2)*N_size) # detector column count (horizontal)
Vert_det = N_size # detector row count (vertical) (no reason for it to be > N)
angles_num = 360 # angles number
angles = np.linspace(0.0,179.9,angles_num,dtype='float32') # in degrees
angles_rad = angles*(np.pi/180.0)
#%%
print ("Building 3D numerical projection data with ASTRA-toolbox")
from tomophantom.supp.astraOP import AstraTools3D

Atools = AstraTools3D(Horiz_det, Vert_det, angles_rad, N_size) # initiate a class object

projData3D_astra = Atools.forwproj(phantom_tm) # numerical projection data

intens_max = 60
sliceSel = 150
#plt.gray()
plt.figure() 
plt.subplot(131)
plt.imshow(projData3D_astra[:,sliceSel,:],vmin=0, vmax=intens_max)
plt.title('2D Projection (astra)')

plt.subplot(132)
plt.imshow(projData3D_astra[sliceSel,:,:],vmin=0, vmax=intens_max)
plt.title('Sinogram view')

plt.subplot(133)
plt.imshow(projData3D_astra[:,:,sliceSel],vmin=0, vmax=intens_max)
plt.title('Tangentogram view')
plt.show()
#%%
# projData3D_astra = np.load('/home/kjy41806/Documents/SOFT/tempProj.npy')

print ("Building 3D analytical projection data with TomoPhantom")
projData3D_analyt= TomoP3D.ModelSino(model, N_size, Horiz_det, Vert_det, angles, path_library3D)

#data rearranging to fit ASTRAs conventions
projData3D_analyt_r = np.zeros((Vert_det, angles_num, Horiz_det),'float32')
for i in range(0,Horiz_det): 
    projData3D_analyt_r[:,:,i] = np.transpose(projData3D_analyt[:,:,i])
projData3D_analyt = projData3D_analyt_r
del projData3D_analyt_r

intens_max = 60
sliceSel = 150
plt.figure() 
plt.subplot(131)
plt.imshow(projData3D_analyt[:,sliceSel,:],vmin=0, vmax=intens_max)
plt.title('2D Projection (analytical)')
plt.subplot(132)
plt.imshow(projData3D_analyt[sliceSel,:,:],vmin=0, vmax=intens_max)
plt.title('Sinogram view')
plt.subplot(133)
plt.imshow(projData3D_analyt[:,:,sliceSel],vmin=0, vmax=intens_max)
plt.title('Tangentogram view')
plt.show()

intens_max = 2
plt.figure() 
plt.subplot(131)
plt.imshow(abs(projData3D_analyt[:,sliceSel,:] - projData3D_astra[:,sliceSel,:]),vmin=0, vmax=intens_max)
plt.title('2D Projection differnce')
plt.subplot(132)
plt.imshow(abs(projData3D_analyt[sliceSel,:,:] - projData3D_astra[sliceSel,:,:]) ,vmin=0, vmax=intens_max)
plt.title('Sinogram difference')
plt.subplot(133)
plt.imshow(abs(projData3D_analyt[:,:,sliceSel] - projData3D_astra[:,:,sliceSel]),vmin=0, vmax=intens_max)
plt.title('Tangentogram difference')
plt.show()
#%%
print ("Reconstruction using ASTRA-toolbox")
recNumerical= Atools.cgls3D(projData3D_analyt, 10) # reconstruct projection data

sliceSel = int(0.5*N_size)
max_val = 1
#plt.gray()
plt.figure(1) 
plt.subplot(131)
plt.imshow(recNumerical[sliceSel,:,:],vmin=0, vmax=max_val)
plt.title('3D Reconstruction, axial view')

plt.subplot(132)
plt.imshow(recNumerical[:,sliceSel,:],vmin=0, vmax=max_val)
plt.title('3D Reconstruction, coronal view')

plt.subplot(133)
plt.imshow(recNumerical[:,:,sliceSel],vmin=0, vmax=max_val)
plt.title('3D Reconstruction, sagittal view')
plt.show()

#%%
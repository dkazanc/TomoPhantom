#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

* Script to generate 3D analytical phantoms and their projection data using TomoPhantom
* Projection data is also generated numerically and reconstructed using ASTRA TOOLBOX 

>>>>> Prerequisites: ASTRA toolbox  <<<<<
install ASTRA: conda install -c astra-toolbox astra-toolbox

@author: Daniil Kazantsev
"""
import timeit
import os
import matplotlib.pyplot as plt
import numpy as np
import tomophantom
from tomophantom import TomoP3D
from tomophantom.supp.qualitymetrics import QualityTools

print ("Building 3D phantom using TomoPhantom software")
tic=timeit.default_timer()
model = 13 # select a model number from the library
N_size = 256 # Define phantom dimensions using a scalar value (cubic phantom)
path = os.path.dirname(tomophantom.__file__)
path_library3D = os.path.join(path, "Phantom3DLibrary.dat")
#This will generate a N_size x N_size x N_size phantom (3D)
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

# Projection geometry related parameters:
Horiz_det = int(np.sqrt(2)*N_size) # detector column count (horizontal)
Vert_det = N_size # detector row count (vertical) (no reason for it to be > N)
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0.0,179.9,angles_num,dtype='float32') # in degrees
angles_rad = angles*(np.pi/180.0)
#%%
print ("Building 3D numerical projection data with ASTRA-toolbox")
from tomophantom.supp.astraOP import AstraTools3D

Atools = AstraTools3D(Horiz_det, Vert_det, angles_rad, N_size) # initiate a class object

projData3D_astra = Atools.forwproj(phantom_tm) # numerical projection data

intens_max = 70
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
print ("Building 3D analytical projection data with TomoPhantom")
projData3D_analyt= TomoP3D.ModelSino(model, N_size, Horiz_det, Vert_det, angles, path_library3D)

intens_max = 70
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

"""
# comparing numerical projections with analytical ones
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
"""
#%%
print ("Reconstruction using ASTRA-toolbox")
recNumerical= Atools.cgls3D(projData3D_analyt, 10) # CGLS-reconstruct projection data

sliceSel = int(0.5*N_size)
max_val = 1
#plt.gray()
plt.figure() 
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

# calculate errors 
Qtools = QualityTools(phantom_tm, recNumerical)
RMSE = Qtools.rmse()
print("Root Mean Square Error is {}".format(RMSE))
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing with FISTA method (ASTRA used for projection)")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

# install FISTA-tomo with: conda install -c dkazanc fista-tomo
# or from https://github.com/dkazanc/FISTA-tomo

from fista.tomo.recModIter import RecTools

# set geometry related parameters and initiate a class object
Rectools = RecTools(DetectorsDimH = Horiz_det,  # DetectorsDimH # detector dimension (horizontal)
                    DetectorsDimV = Vert_det,  # DetectorsDimV # detector dimension (vertical) for 3D case only
                    AnglesVec = angles_rad, # array of angles in radians
                    ObjSize = N_size, # a scalar to define reconstructed object dimensions
                    datafidelity='LS',# data fidelity, choose LS, PWLS (wip), GH (wip), Student (wip)
                    OS_number = None, # the number of subsets, NONE/(or > 1) ~ classical / ordered subsets
                    tolerance = 1e-06, # tolerance to stop outer iterations earlier
                    device='gpu')


lc = Rectools.powermethod() # calculate Lipschitz constant

# Run FISTA reconstrucion algorithm without regularisation
RecFISTA = Rectools.FISTA(projData3D_analyt, iterationsFISTA = 85, lipschitz_const = lc)

# Run FISTA reconstrucion algorithm with 3D regularisation
#RecFISTA_reg = Rectools.FISTA(projData3D_analyt, iterationsFISTA = 85, regularisation = 'ROF_TV', lipschitz_const = lc)

sliceSel = int(0.5*N_size)
max_val = 1
plt.figure() 
plt.subplot(131)
plt.imshow(RecFISTA[sliceSel,:,:],vmin=0, vmax=max_val)
plt.title('3D FISTA Reconstruction, axial view')

plt.subplot(132)
plt.imshow(RecFISTA[:,sliceSel,:],vmin=0, vmax=max_val)
plt.title('3D FISTA Reconstruction, coronal view')

plt.subplot(133)
plt.imshow(RecFISTA[:,:,sliceSel],vmin=0, vmax=max_val)
plt.title('3D FISTA Reconstruction, sagittal view')
plt.show()
#%%



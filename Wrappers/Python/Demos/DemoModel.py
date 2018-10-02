#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D analytical phantoms and their sinograms
If one needs to modify/add phantoms, please edit Phantom2DLibrary.dat

>>>>> Optional prerequisites: ASTRA toolbox <<<<<
install ASTRA: conda install -c astra-toolbox astra-toolbox
 
!Run this script from "Demos" folder in order to ensure a correct path to *dat file!

@author: Daniil Kazantsev
"""
import numpy as np
import timeit
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
plt.figure()
plt.rcParams.update({'font.size': 21})
plt.imshow(phantom_2D, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}''{}'.format('2D Phantom using model no.',model))

# parameters to generate a sinogram
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0.0,179.9,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180.0)
P = int(np.sqrt(2)*N_size) #detectors
#%%
###################################################################
tic=timeit.default_timer()
# create sinogram analytically
sino_an = TomoP2D.ModelSino(model, N_size, P, angles, pathTP)
toc=timeit.default_timer()
Run_time = toc - tic
print("Analytical sinogram has been generated in {} seconds".format(Run_time))

plt.figure()
plt.rcParams.update({'font.size': 21})
plt.imshow(sino_an,  cmap="BuPu")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Analytical sinogram of model no.',model))
#%%
###################################################################
# get numerical sinogram
tic=timeit.default_timer()
sino_num = TomoP2D.SinoNum (phantom_2D, P, angles)
toc=timeit.default_timer()
Run_time = toc - tic
print("Numerical sinogram has been generated in {} seconds".format(Run_time))

plt.figure()
plt.rcParams.update({'font.size': 21})
plt.imshow(sino_num,  cmap="BuPu")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Numerical sinogram of model no.',model))
#%%
###################################################################
# get numerical sinogram (ASTRA-toolbox)
from tomophantom.supp.astraOP import AstraTools

tic=timeit.default_timer()
Atools = AstraTools(P, angles_rad, N_size, 'cpu') # initiate a class object
sino_num_ASTRA = Atools.forwproj(phantom_2D) # generate numerical sino (Ax)
toc=timeit.default_timer()
Run_time = toc - tic
print("Numerical (ASTRA) sinogram has been generated in {} seconds".format(Run_time))

plt.figure()
plt.rcParams.update({'font.size': 21})
plt.imshow(sino_num_ASTRA,  cmap="BuPu")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Numerical sinogram (ASTRA) of model no.',model))
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using Fourier Slice method")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
from tomophantom.supp.recMod import RecTools

Rectools = RecTools(P, angles_rad, N_size) # initiate a class object
RecFourier = Rectools.fourier(sino_an,'linear') 

plt.figure() 
plt.imshow(RecFourier, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('Fourier slice reconstruction')
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using FBP (ASTRA-TOOLBOX)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
from tomophantom.supp.astraOP import AstraTools

#x = Atools.backproj(sino_an) # generate backprojection (A'b)

plt.figure() 
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
plt.figure() 
plt.imshow(FBPrec1, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP Reconstructed Phantom (analyt)')

print ("Reconstructing numerical sinogram using FBP (astra TB)...")
FBPrec2 = Atools.fbp2D(sino_num_ASTRA)

plt.figure() 
plt.imshow(FBPrec2, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP Reconstructed Phantom (numeric)')

plt.figure() 
plt.imshow(abs(FBPrec1-FBPrec2), vmin=0, vmax=0.05, cmap="BuPu")
plt.colorbar(ticks=[0, 0.02, 0.05], orientation='vertical')
plt.title('FBP rec differences')
rmse2 = np.linalg.norm(FBPrec1 - FBPrec2)/np.linalg.norm(FBPrec2)

#%%

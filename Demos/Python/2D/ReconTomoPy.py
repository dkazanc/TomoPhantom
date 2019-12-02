#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D analytical phantoms and their sinograms with added noise and artifacts
Sinograms then reconstructed using TomoPy package

>>>> Prerequisites: TomoPy package <<<<<

This demo demonstrates frequent inaccuracies which are accosiated with X-ray imaging:
zingers, rings and noise

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP2D
import os
import tomophantom

model = 4 # select a model
N_size = 512
# one can specify an exact path to the parameters file
# path_library2D = '../../../PhantomLibrary/models/Phantom2DLibrary.dat'
path = os.path.dirname(tomophantom.__file__)
path_library2D = os.path.join(path, "Phantom2DLibrary.dat")
phantom_2D = TomoP2D.Model(model, N_size, path_library2D)

plt.close('all')
plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.imshow(phantom_2D, vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}''{}'.format('2D Phantom using model no.',model))

# create sinogram analytically
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors

sino_an = TomoP2D.ModelSino(model, N_size, P, angles, path_library2D)

plt.figure(2)
plt.rcParams.update({'font.size': 21})
plt.imshow(sino_an,  cmap="gray")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Analytical sinogram of model no.',model))

#%%
# Adding noise
from tomophantom.supp.artifacts import _Artifacts_

# forming dictionaries with artifact types
_noise_ =  {'type' : 'Poisson',
            'sigma' : 10000, # noise amplitude
            'seed' : 0}

noisy_sino = _Artifacts_(sino_an, _noise_, _zingers_={}, _stripes_={}, _sinoshifts_= {})

plt.figure()
plt.rcParams.update({'font.size': 21})
plt.imshow(noisy_sino,cmap="gray")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}''{}'.format('Analytical noisy sinogram.',model))
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using gridrec (TomoPy)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
import tomopy

sinoTP = np.zeros((angles_num,1,P),dtype='float32')
sinoTP[:,0,:] = noisy_sino
rot_center = tomopy.find_center(sinoTP, angles_rad, init=290, ind=0, tol=0.5)
reconTomoPy_ideal = tomopy.recon(sinoTP, angles_rad, center=rot_center, algorithm='gridrec')

plt.figure()
plt.imshow(reconTomoPy_ideal[0,:,:], vmin=0, vmax=1, cmap="gray")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('GridRec reconstruction (TomoPy)')
plt.show()
#%%
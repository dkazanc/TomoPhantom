#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 20:38:35 2018

@author: algol
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D/3D analytical phantoms and their sinograms
If one needs to modify/add phantoms, please edit Phantom2DLibrary.dat or
Phantom3DLibrary.dat
>>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

@author: daniil.kazantsev@manchester.ac.uk
"""
import astra
import numpy as np
import matplotlib.pyplot as plt
import sys 
import os

# modify your path to the compiled package accordingly
#sys.path.append('~/.local/lib/python3.6/site-packages/')
path = os.getcwd()
head, tail = os.path.split(path)
path = head
sys.path.append(path)
from tomophantom import phantom2d
from tomophantom import phantom3d

#%%
model = 11
N_size = 512
pathTP = '../../functions/models/Phantom2DLibrary.dat'
#This will generate a N_size x N_size phantom (2D)
phantom_2D = phantom2d.buildPhantom2D(model, N_size, pathTP)

plt.figure(1)
plt.imshow(phantom_2D, cmap="hot")
plt.title('2D Phantom')
#%%
angles_num = int(np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors

# create sinogram analytically
sino_an = phantom2d.buildSino2D(model, N_size, P, angles, pathTP, 1)

plt.figure(2)
plt.imshow(sino_an, cmap="hot")
plt.title('Analytical sinogram')
#%%
# create sinogram using ASTRA toolbox
proj_geom = astra.create_proj_geom('parallel', 1.0, P, angles_rad - 0.5*np.pi)
vol_geom = astra.create_vol_geom(N_size, N_size)
# Create 2D projection data
proj_id = astra.create_projector('cuda',proj_geom,vol_geom)
sinogram_id, sino_num_ASTRA = astra.create_sino(phantom_2D, proj_id)
astra.data2d.delete(sinogram_id)

plt.figure(2) 
plt.subplot(121)
plt.imshow(sino_an, cmap="hot")
plt.title('Analytical sinogram')
plt.subplot(122)
plt.imshow(sino_num_ASTRA, cmap="hot")
plt.title('Numerical sinogram')
plt.show()
#%%
print ("Reconstructing...")

rec_id = astra.data2d.create( '-vol', vol_geom)

# Create a data object to hold the sinogram data
sinogram_id = astra.data2d.create('-sino', proj_geom, sino_an)

cfg = astra.astra_dict('FBP_CUDA')
cfg['ReconstructionDataId'] = rec_id
cfg['ProjectionDataId'] = sinogram_id
cfg['FilterType'] = 'Ram-Lak'

# Create and run the algorithm object from the configuration structure
alg_id = astra.algorithm.create(cfg)
astra.algorithm.run(alg_id)
# Get the result
rec = astra.data2d.get(rec_id)

astra.algorithm.delete(alg_id)
astra.data2d.delete(rec_id)
astra.data2d.delete(sinogram_id)

plt.figure(3) 
plt.imshow(rec, vmin=0, vmax=1, cmap="hot")
plt.title('Reconstructed Phantom')

#%%
print ("Building 3D phantom using TomoPhantom software")
model = 7
N_size = 512
pathTP3 = '../../functions/models/Phantom3DLibrary.dat'
#This will generate a N_size x N_size x N_size phantom (3D)
phantom_tm = phantom3d.buildPhantom3D(model, N_size, pathTP3)

#plt.gray()
plt.figure(4) 
plt.subplot(121)
plt.imshow(phantom_tm[128,:,:],vmin=0, vmax=1, cmap="hot")
plt.title('3D Phantom, axial view')

plt.subplot(122)
plt.imshow(phantom_tm[:,128,:],vmin=0, vmax=1, cmap="hot")
plt.title('3D Phantom, coronal view')
plt.show()


"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate temporal (3D and 4D) analytical phantoms and their sinograms
If one needs to modify/add phantoms, please edit Phantom2DLibrary.dat or
Phantom3DLibrary.dat
Note that all temporal phantoms start from no. 100
>>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

Run demo from the folder "Demos"

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP2D

model = 102  # note that the selected model is temporal (2D + time)
N_size = 512
timeframes = 25
#specify a full path to the parameters file
pathTP = '../../functions/models/Phantom2DLibrary.dat'
#This will generate a N_size x N_size x Time frames phantom (2D + time)
phantom_2Dt = TomoP2D.ModelTemporal(model, N_size, pathTP)

plt.close('all')
plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.title('{}''{}'.format('2D+t phantom using model no.',model))
for sl in range(0,timeframes):
    im = phantom_2Dt[sl,:,:]
    plt.imshow(im, vmin=0, vmax=1)
    plt.pause(.1)
    plt.draw

# create sinogram analytically
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors

sino = TomoP2D.ModelSinoTemporal(model, N_size, P, angles, pathTP)

plt.figure(2)
plt.rcParams.update({'font.size': 21})
plt.title('{}''{}'.format('2D+t sinogram of model no.',model))
for sl in range(0,timeframes):
    im = sino[sl,:,:].transpose()
    plt.imshow(im, vmin=0, vmax=180)
    plt.pause(.1)
    plt.draw
#%%
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print ("Reconstructing analytical sinogram using FBP (ASTRA-TOOLBOX)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
from tomophantom.supp.astraOP import AstraTools
Atools = AstraTools(P, angles_rad, N_size, 'cpu') # initiate a class object
FBPrec = Atools.fbp2D(sino[15,:,:].transpose())

plt.figure(3) 
plt.imshow(FBPrec, vmin=0, vmax=1)
plt.title('FBP Reconstructed Phantom')
#%%
"""
from tomophantom import TomoP3D

# generate 4D (3D + time) model
model = 101 # note that the selected model is temporal (3D + time)
N_size = 256
#specify a full path to the parameters file
pathTP = '../../functions/models/Phantom3DLibrary.dat'
#This will generate a Time frames x N_size x N_size N_size phantom (3D + time)
phantom_3Dt = TomoP3D.ModelTemporal(model, N_size, pathTP)


sliceSel = int(0.5*N_size)
#plt.gray()
plt.figure(4) 
plt.imshow(phantom_3Dt[0,sliceSel,:,:],vmin=0, vmax=1)
plt.title('4D Phantom, axial view, first time-frame')
"""
#%%
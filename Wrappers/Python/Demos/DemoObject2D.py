#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D/3D analytical objects and their sinograms
Recursively adding objects and sinos one can build a required model

>>>>> Optional dependencies (reconstruction mainly): <<<<<
1. ASTRA toolbox: conda install -c astra-toolbox astra-toolbox
2. TomoRec: conda install -c dkazanc tomorec
or install from https://github.com/dkazanc/TomoRec

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP2D 
from tomophantom.TomoP2D import Objects2D

# create a 2D object explicitly without using parameters file
N_size = 512 # define the grid

# define all objects bellow:
pp = {'Obj': Objects2D.GAUSSIAN,
      'C0' : 1.00,
      'x0' : 0.25,
      'y0' : -0.3,
      'a'  : 0.15,
      'b'  :  0.3,
      'phi': -30.0}

pp1 = {'Obj': Objects2D.RECTANGLE,
      'C0' : 1.00,
      'x0' : -0.2,
      'y0' : 0.2,
      'a'  : 0.25,
      'b'  :  0.4,
      'phi': 60.0}

myObjects = [pp, pp1] # dictionary of objects
Object1 = TomoP2D.Object(N_size, myObjects)

plt.close('all')
plt.figure(1)
plt.rcParams.update({'font.size': 21})
plt.imshow(Object1, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('{}'.format('2D Object'))

# create sinogram analytically without using the parameters file
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors

sino_an = TomoP2D.ObjectSino(N_size, P, angles, myObjects)

plt.figure(2)
plt.rcParams.update({'font.size': 21})
plt.imshow(sino_an, cmap="BuPu")
plt.colorbar(ticks=[0, 150, 250], orientation='vertical')
plt.title('{}'.format('Analytical sinogram of an object'))
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

FBPrec = RectoolsDIR.FBP(sino_an)

plt.figure(3) 
plt.imshow(FBPrec, vmin=0, vmax=1, cmap="BuPu")
plt.title('FBP Reconstructed Model')
#%%
"""
Main difference from DemoModel.py is that we extract all parameters from the 
library file using Python and then pass it to the Object function instead of model.
This can be helpful if one would like to avoid using library files and can
pass parameters directly into object function

!Run script from "Demos" folder in order to ensure a correct path to *dat file!
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import tomophantom
from tomophantom import TomoP2D
from tomophantom.supp.libraryToDict import modelfile2Dtolist

model = 11 # select a model number from the library
N_size = 512 # set dimension of the phantom
# one can specify an exact path to the parameters file
# path_library2D = '../../../PhantomLibrary/models/Phantom2DLibrary.dat'
path = os.path.dirname(tomophantom.__file__)
path_library2D = os.path.join(path, "Phantom2DLibrary.dat")
objlist = modelfile2Dtolist(path_library2D, model) # extract parameters using Python
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
print ("Reconstructing analytical sinogram using FBP (TomoRec)...")
print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
# initialise TomoRec reconstruction class ONCE
from tomorec.methodsDIR import RecToolsDIR
RectoolsDIR = RecToolsDIR(DetectorsDimH = P,  # DetectorsDimH # detector dimension (horizontal)
                    DetectorsDimV = None,  # DetectorsDimV # detector dimension (vertical) for 3D case only
                    AnglesVec = angles_rad, # array of angles in radians
                    ObjSize = N_size, # a scalar to define reconstructed object dimensions
                    device='cpu')


sino_num_ASTRA = RectoolsDIR.FORWPROJ(phantom_2D) # generate numerical sino (Ax)
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
FBPrec1 = RectoolsDIR.FBP(sino_an)
plt.figure(4) 
plt.imshow(FBPrec1, vmin=0, vmax=1, cmap="BuPu")
plt.colorbar(ticks=[0, 0.5, 1], orientation='vertical')
plt.title('FBP Reconstructed Phantom (analyt)')

print ("Reconstructing numerical sinogram using FBP (astra TB)...")
FBPrec2 = RectoolsDIR.FBP(sino_num_ASTRA)
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


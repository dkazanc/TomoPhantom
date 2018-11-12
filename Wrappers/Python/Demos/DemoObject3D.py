#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 3D analytical objects and their projection data
Recursively adding objects one can build a required model with the corresponding projection data

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP3D
from tomophantom.TomoP3D import Objects3D

N3D_size = 256

# specify object parameters, here we replicate model 
obj3D_1 = {'Obj': Objects3D.GAUSSIAN, 
      'C0' : 1.00,
      'x0' :-0.25,
      'y0' : -0.15,
      'z0' : 0.0,
      'a'  : 0.3,
      'b'  :  0.2,
      'c'  :  0.3,
      'phi1'  : 35.0}

obj3D_2 = {'Obj': Objects3D.CUBOID, 
      'C0' : 1.00,
      'x0' :0.1,
      'y0' : 0.2,
      'z0' : 0.0,
      'a'  : 0.15,
      'b'  :  0.35,
      'c'  :  0.6,
      'phi1'  : -60.0}

print ("Building 3D object using TomoPhantom software")
myObjects = [obj3D_1, obj3D_2] # dictionary of objects
Object3D = TomoP3D.Object(N3D_size, myObjects)

sliceSel = int(0.5*N3D_size)
#plt.gray()
plt.figure() 
plt.subplot(131)
plt.imshow(Object3D[sliceSel,:,:],vmin=0, vmax=1)
plt.title('3D Object, axial view')

plt.subplot(132)
plt.imshow(Object3D[:,sliceSel,:],vmin=0, vmax=1)
plt.title('3D Object, coronal view')

plt.subplot(133)
plt.imshow(Object3D[:,:,sliceSel],vmin=0, vmax=1)
plt.title('3D Object, sagittal view')
plt.show()
#%%
Horiz_det = int(np.sqrt(2)*N3D_size) # detector column count (horizontal)
Vert_det = N3D_size # detector row count (vertical) (no reason for it to be > N)
angles_num = int(0.5*np.pi*N3D_size); # angles number
angles = np.linspace(0.0,179.9,angles_num,dtype='float32') # in degrees

print ("Building 3D analytical projection data with TomoPhantom")
ProjData3D = TomoP3D.ObjectSino(N3D_size, Horiz_det, Vert_det, angles, myObjects)

intens_max = 60
sliceSel = 150
plt.figure() 
plt.subplot(131)
plt.imshow(ProjData3D[:,sliceSel,:],vmin=0, vmax=intens_max)
plt.title('2D Projection (analytical)')
plt.subplot(132)
plt.imshow(ProjData3D[sliceSel,:,:],vmin=0, vmax=intens_max)
plt.title('Sinogram view')
plt.subplot(133)
plt.imshow(ProjData3D[:,:,sliceSel],vmin=0, vmax=intens_max)
plt.title('Tangentogram view')
plt.show()
#%%



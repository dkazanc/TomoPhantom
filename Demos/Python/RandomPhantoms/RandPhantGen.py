#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The TomoPhantom package is released under Apache License, Version 2.0

Script to generate 2D/3D random objects and their projection data
- foam2D/3D generate a phantom with non-overlapping spherical objects,
for 2D object_type choose: 'ellipse', 'parabola', 'gaussian' or 'mix'
and for 3D: 'ellipsoid', 'paraboloid', 'gaussian', 'mix'

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt

from tomophantom import TomoP2D 
from tomophantom import TomoP3D 
from tomophantom.randphant.generator import foam2D,foam3D

N_size = 256 # define the grid size
tot_objects = 300 # the total number of objects to generate

# define ranges for parameters
x0min = -0.9
x0max = 0.9
y0min = -0.9
y0max = 0.9
z0min = -0.9
z0max = 0.9
c0min = 0.01
c0max = 1.0
ab_min = 0.01
ab_max = 0.25

# 2D example
(Objfoam2D,myObjects) = foam2D(x0min, x0max, y0min, y0max, c0min, c0max, ab_min, ab_max, N_size, tot_objects, object_type = 'mix')
plt.figure()
plt.imshow(Objfoam2D,  vmin=0, vmax=0.3, cmap="gray")

# Generate a sinogram
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0,180,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180)
P = int(np.sqrt(2)*N_size) #detectors

sino_Objfoam2D = TomoP2D.ObjectSino(N_size, P, angles, myObjects)
plt.figure()
plt.imshow(sino_Objfoam2D,  cmap="gray")

#%%
# 3D example
print ("Generating 3D random phantom")
(Objfoam3D, myObjects) = foam3D(x0min, x0max, y0min, y0max, z0min, z0max, c0min, c0max, ab_min, ab_max, N_size, tot_objects, object_type = 'mix')
plt.figure()
sliceSel = int(0.5*N_size)
plt.subplot(131)
plt.imshow(Objfoam3D[sliceSel,:,:],vmin=0, vmax=0.5)
plt.title('3D Object, axial view')

plt.subplot(132)
plt.imshow(Objfoam3D[:,sliceSel,:],vmin=0, vmax=0.5)
plt.title('3D Object, coronal view')

plt.subplot(133)
plt.imshow(Objfoam3D[:,:,sliceSel],vmin=0, vmax=0.5)
plt.title('3D Object, sagittal view')
plt.show()

Horiz_det = int(np.sqrt(2)*N_size) # detector column count (horizontal)
Vert_det = N_size # detector row count (vertical) (no reason for it to be > N)
angles_num = int(0.5*np.pi*N_size); # angles number
angles = np.linspace(0.0,179.9,angles_num,dtype='float32') # in degrees

print ("Building 3D analytical projection data with TomoPhantom")
ProjData3D = TomoP3D.ObjectSino(N_size, Horiz_det, Vert_det, angles, myObjects)

sliceSel = 150
plt.figure() 
plt.subplot(131)
plt.imshow(ProjData3D[:,sliceSel,:])
plt.title('2D Projection (analytical)')
plt.subplot(132)
plt.imshow(ProjData3D[sliceSel,:,:])
plt.title('Sinogram view')
plt.subplot(133)
plt.imshow(ProjData3D[:,:,sliceSel])
plt.title('Tangentogram view')
plt.show()
#%%
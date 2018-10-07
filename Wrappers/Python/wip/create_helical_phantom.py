#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 15:21:35 2018
Script to generate helical scan data

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP3D
from tomophantom.TomoP3D import Objects3D
from tomophantom.supp.astraOP import AstraTools3D

angles_num = 20 # the total number of projection angles
angles_range_min = 0 # acquisition range (start)
angles_range_max = 720 # acquisition range (finish)
angles = np.linspace(angles_range_min,angles_range_max,angles_num,dtype='float32')
angles_rad = angles*(np.pi/180.0) # obtained radian angles

scale_factor = 10 # factor to scale the problem down 

# detector dimensions
x_dim = int(2000/scale_factor)
y_dim = int(2000/scale_factor)
z_dim = int(2540/scale_factor)
z_step = 0.7/scale_factor # step to increase vertical step of objects
z_start = int(450/scale_factor) # initial positioning

steps_no = round((z_dim - z_start)/z_step) # the total number of steps
#converting values to the coordinate system of an object
z0 = z_start/z_dim
z_start_obj = -(1.0 - 2.0*z0) # starting z coordinate of objects
z_coord_steps = np.linspace(z_start_obj,1.0,steps_no,dtype='float32')

ObjSize = (z_dim,y_dim,x_dim) # set dimensions of the object
DetRows = ObjSize[0] # detector vertical dimension (rows)
DetColumns = ObjSize[1] # detector horizontal dimension (columns)
proj3D = np.zeros((angles_num,DetRows,DetColumns),'float32')

# loop over angles
for i in range(0,angles_num):
    obj3D1 = {'Obj': Objects3D.ELLIPSOID, 
          'C0' : 1.00,
          'x0' :-0.5,
          'y0' : 0.5,
          'z0' : z_coord_steps[i],
          'a'  : 0.15,
          'b'  :  0.15,
          'c'  :  0.15,
          'phi1'  :  0.0,
          'phi2'  :  0.0,
          'phi3': 0.0}
    obj3D2 = {'Obj': Objects3D.ELLIPSOID, 
          'C0' : 1.00,
          'x0' : 0.5,
          'y0' : 0.5,
          'z0' : z_coord_steps[i],
          'a'  : 0.15,
          'b'  :  0.1,
          'c'  :  0.2,
          'phi1'  :  0.0,
          'phi2'  :  0.0,
          'phi3': 0.0}
    obj3D3 = {'Obj': Objects3D.ELLIPSOID, 
          'C0' : 1.00,
          'x0' : 0.5,
          'y0' : -0.5,
          'z0' : z_coord_steps[i],
          'a'  : 0.2,
          'b'  :  0.15,
          'c'  :  0.1,
          'phi1'  :  0.0,
          'phi2'  :  0.0,
          'phi3': 0.0}
    obj3D4 = {'Obj': Objects3D.ELLIPSOID, 
          'C0' : 1.00,
          'x0' : -0.5,
          'y0' : -0.5,
          'z0' : z_coord_steps[i],
          'a'  : 0.18,
          'b'  :  0.07,
          'c'  :  0.13,
          'phi1'  :  0.0,
          'phi2'  :  0.0,
          'phi3': 0.0}
        
    myObjects = [obj3D1, obj3D2, obj3D3, obj3D4] # dictionary of objects
    
    object3D = TomoP3D.Object(ObjSize, myObjects)
    """
    sliceSel = 200
    #plt.gray()
    plt.figure() 
    plt.subplot(121)
    plt.imshow(Object3D[sliceSel,:,:],vmin=0, vmax=1)
    plt.title('3D Phantom, axial view')
    
    plt.subplot(122)
    plt.imshow(Object3D[:,sliceSel,:],vmin=0, vmax=1)
    plt.title('3D Phantom, coronal view')
    plt.show()
    """
    # obtain projections using ASTRA-toolbox
    Atools = AstraTools3D(DetRows, DetColumns, angles_rad, ObjSize) # initiate a class object
    projection2D = Atools.forwproj(object3D) # get a projection for specific angle
    
    #collecting into 3D projection data from each angle
    proj3D[i,:,:] = projection2D[:,0,:]

"""
for i in range(0,angles_num):
    plt.figure(10) 
    plt.imshow(proj3D[i,:,:],vmin=0, vmax=100)
    plt.show()
    plt.pause(0.1)
"""
#%%
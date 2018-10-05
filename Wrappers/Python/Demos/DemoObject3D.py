#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPLv3 license (ASTRA toolbox)
Note that the TomoPhantom package is released under Apache License, Version 2.0

Script to generate 3D analytical objects and their sinograms
Recursively adding objects and sinos one can build a required model

>>>> Prerequisites: ASTRA toolbox, if one needs to do reconstruction <<<<<

Run demo from the folder "Demos"

@author: Daniil Kazantsev
"""
import numpy as np
import matplotlib.pyplot as plt
from tomophantom import TomoP3D
from tomophantom.TomoP3D import Objects3D

N3D_size = 256

obj3D = {'Obj': Objects3D.GAUSSIAN, 
      'C0' : 1.00,
      'x0' :-0.25,
      'y0' : 0.1,
      'z0' : 0.0,
      'a'  : 0.2,
      'b'  :  0.35,
      'c'  :  0.7,
      'phi1'  :  30.0,
      'phi2'  :  60.0,
      'phi3': -25.0}

Object3D = TomoP3D.Object(N3D_size, [obj3D])

sliceSel = int(0.5*N3D_size)
#plt.gray()
plt.figure(7) 
plt.subplot(121)
plt.imshow(Object3D[sliceSel,:,:],vmin=0, vmax=1)
plt.title('3D Phantom, axial view')

plt.subplot(122)
plt.imshow(Object3D[:,sliceSel,:],vmin=0, vmax=1)
plt.title('3D Phantom, coronal view')
plt.show()
#%%


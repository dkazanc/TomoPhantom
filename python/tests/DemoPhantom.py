#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 10:45:42 2017

@author: algol
"""
#import numpy as np
import pylab

from tomophantom import phantom3d

model = 7
N_size = 512
path = '/home/algol/Documents/MATLAB/TomoPhantom/functions/models/Phantom3DLibrary.dat'
#This will generate a N_size x N_size x N_size phantom (3D)
phantom_tm = phantom3d.build_volume_phantom_3d(path, model, N_size)

pylab.figure(1) 
pylab.subplot(211)
pylab.imshow(phantom_tm[300,:,:],vmin=0, vmax=1)
pylab.title('3D Phantom, axial view')

pylab.subplot(212)
pylab.imshow(phantom_tm[:,300,:],vmin=0, vmax=1)
pylab.title('3D Phantom, coronal view')
pylab.show()  
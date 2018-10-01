# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 16:21:58 2018

@author: ofn77899
"""

from tomophantom import TomoP3D as tp3
import matplotlib.pyplot as plt
import numpy as np

#tp3.Object(tp3.)

# size of phantom
n1=n2=n3 = 128
size = (n1,n2,n3)

# number of objects
M = 10

#margin
margin = 3

# random location
X = list(np.random.uniform(low = 0, high=1, size=M))
Y = np.asarray(np.random.uniform(low = 0, high=1, size=M), dtype=float)
Z = np.asarray(np.random.uniform(low = 0, high=1, size=M), dtype=float)




gaussians = [{'Obj': tp3.Objects3D.GAUSSIAN, 
            'C0' : 1., 
            'x0' : X[loc], 
            'y0' : Y[loc], 
            'z0' : 0.,
            'a'  : .2, 
            'b'  : .2 , 
            'c'  : .2,
            'phi1':0., 
            'phi2':0.,
            'phi3':0. } for loc in range(len(X))]

vol = tp3.Object(128, gaussians[0])
plt.imshow(vol[0])
plt.show()

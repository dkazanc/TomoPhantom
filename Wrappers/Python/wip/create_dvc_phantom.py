# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 16:21:58 2018

@author: ofn77899
"""

from tomophantom import TomoP3D as tp3
import matplotlib.pyplot as plt
import numpy as np
import vtk
from ccpi.viewer.CILViewer2D import Converter

import time

#tp3.Object(tp3.)

class NonOverlappingPoints(object):
    def __init__(self, minimum_distance):
        self.minimum_distance = minimum_distance
        self.points = []
        self.apoints = None
        self.rejected = 0
    
    @staticmethod
    def distance(point1 , point2):
        vec = point1 - point2
        return np.sqrt(np.dot(vec.T, vec))
    
    def append(self, point):
        if self.apoints is None:
            self.apoints = np.asarray(point, dtype=np.float32)
            return 1
        
        if self.apoints.size > 3:
            b = self.apoints - point
            diag = np.diagonal(np.dot(b,b.T))
            m = np.sort(diag)[-1]
            N = self.apoints.shape[0]
        else:
            m = np.dot(self.apoints, point)
            N = 1
        if m < self.minimum_distance*self.minimum_distance:
            self.rejected += 1
            return N
        else:
            self.apoints = np.vstack((self.apoints, point))
            return self.apoints.shape[0]
    
    def append_for(self, point):
        M = len (self.points)
        if M == 0:
            self.points.append(np.asarray(point, dtype=np.float32))
            return len(self.points)
        overlaps = False
        for i in range (M):
            if NonOverlappingPoints.distance(point, self.points[i]) < \
               self.minimum_distance:
                   overlaps = True
                   self.rejected += 1
                   break
           
        if not overlaps:
            self.points.append(np.asarray(point, dtype=np.float32))
        
        return len (self.points)
        



# size of phantom
n1=n2=n3 = 128
size = (n1,n2,n3)

# number of objects
M = 300

#margin
margin = 3


# random location
X = np.asarray(np.random.uniform(low = -1, high=1, size=M), dtype=np.float32)
Y = np.asarray(np.random.uniform(low = -1, high=1, size=M), dtype=np.float32)
Z = np.asarray(np.random.uniform(low = -1, high=1, size=M), dtype=np.float32)

#points = np.vstack((X,Y,Z)).T

#Gaussian FWHM 2.355 sigma
nop = NonOverlappingPoints(0.2)

#nop.append(np.random.uniform(low = 0, high=1, size=3))
counter = 0
start = time.time()
while (nop.append(np.random.uniform(low = -1, high=1, size=3)) < M):
    pass
    #print ("Adding point: " , counter)
    #counter += 1
end1 = time.time()
print ("append", end1-start)
while (nop.append_for(np.random.uniform(low = -1, high=1, size=3)) < M):
    pass
end2 = time.time()
print ("append_for", end2-end1)


#%%
gaussians = [{'Obj': tp3.Objects3D.PARABOLOID, 
            'C0' : 1., 
            'x0' : float(nop.apoints[loc][0]),
            'y0' : float(nop.apoints[loc][1]), 
            'z0' : float(nop.apoints[loc][2]),
            'a'  : .2, 
            'b'  : .1 , 
            'c'  : .15,
            'phi1':np.random.uniform(low = 0, high=360), 
            'phi2':np.random.uniform(low = 0, high=360),
            'phi3':np.random.uniform(low = 0, high=360) } for loc in range(len(nop.apoints))]

vol = tp3.Object(128, gaussians)
plt.imshow(vol[0])
plt.show()


writer = vtk.vtkMetaImageWriter()
conv = Converter.numpy2vtkImporter(vol)
conv.Update()
writer.SetInputData(conv.GetOutput())
writer.SetFileName("dvc_phantom0.mha")
writer.Write()